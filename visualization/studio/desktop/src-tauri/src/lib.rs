use serde_json::Value;
use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::net::{Ipv4Addr, SocketAddr, SocketAddrV4, TcpListener, TcpStream};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::Mutex;
use std::time::{Duration, Instant};

#[derive(Default)]
struct ViewerProcesses(Mutex<HashMap<String, Child>>);

impl Drop for ViewerProcesses {
    fn drop(&mut self) {
        if let Ok(children) = self.0.get_mut() {
            for (run_directory, child) in children.iter_mut() {
                let _ = child.kill();
                let _ = child.wait();
                let _ = fs::remove_file(
                    PathBuf::from(run_directory)
                        .join("control")
                        .join("viewer.attached"),
                );
            }
            children.clear();
        }
    }
}

fn atomic_write(path: &Path, bytes: &[u8]) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    let temporary = PathBuf::from(format!("{}.tmp", path.display()));
    let mut file = fs::File::create(&temporary).map_err(|e| e.to_string())?;
    file.write_all(bytes).map_err(|e| e.to_string())?;
    file.sync_all().map_err(|e| e.to_string())?;
    fs::rename(temporary, path).map_err(|e| e.to_string())
}

#[tauri::command]
fn choose_yaml_file() -> Option<String> {
    rfd::FileDialog::new()
        .add_filter("YAML", &["yaml", "yml"])
        .pick_file()
        .map(|p| p.to_string_lossy().into_owned())
}

#[tauri::command]
fn choose_run_directory() -> Option<String> {
    rfd::FileDialog::new()
        .pick_folder()
        .map(|p| p.to_string_lossy().into_owned())
}

#[tauri::command]
fn load_yaml(path: String) -> Result<Value, String> {
    let text = fs::read_to_string(path).map_err(|e| e.to_string())?;
    let yaml: serde_yaml::Value = serde_yaml::from_str(&text).map_err(|e| e.to_string())?;
    serde_json::to_value(yaml).map_err(|e| e.to_string())
}

#[tauri::command]
fn save_yaml(path: String, value: Value) -> Result<(), String> {
    let yaml: serde_yaml::Value = serde_json::from_value(value).map_err(|e| e.to_string())?;
    let mut text = serde_yaml::to_string(&yaml).map_err(|e| e.to_string())?;
    if !text.ends_with('\n') {
        text.push('\n');
    }
    atomic_write(Path::new(&path), text.as_bytes())
}

#[tauri::command]
fn validate_config(binary: String, config: String) -> Result<String, String> {
    let output = Command::new(binary)
        .args(["--config", &config, "--dry-run"])
        .output()
        .map_err(|e| e.to_string())?;
    if output.status.success() {
        String::from_utf8(output.stdout).map_err(|e| e.to_string())
    } else {
        Err(String::from_utf8_lossy(&output.stderr).trim().to_string())
    }
}

#[tauri::command]
fn start_simulation(
    launcher: String,
    binary: String,
    config: String,
    run_directory: String,
) -> Result<String, String> {
    let run = PathBuf::from(&run_directory);
    let name = run
        .file_name()
        .and_then(|v| v.to_str())
        .ok_or_else(|| "run directory has no valid name".to_string())?;
    let parent = run
        .parent()
        .ok_or_else(|| "run directory has no parent".to_string())?;
    let prefix = parent.join(".atcg3d-control").join(name).join("runner");
    let output = Command::new(launcher)
        .arg(binary)
        .arg(config)
        .arg(prefix)
        .output()
        .map_err(|e| e.to_string())?;
    if output.status.success() {
        Ok(String::from_utf8_lossy(&output.stdout).trim().to_string())
    } else {
        Err(String::from_utf8_lossy(&output.stderr).trim().to_string())
    }
}

#[tauri::command]
fn send_control(run_directory: String, action: String) -> Result<(), String> {
    let file = match action.as_str() {
        "pause" => "pause.request",
        "resume" => "resume.request",
        "checkpoint_stop" => "checkpoint_stop.request",
        "terminate" => "terminate.request",
        _ => return Err("unsupported control action".to_string()),
    };
    let path = PathBuf::from(run_directory).join("control").join(file);
    atomic_write(&path, b"1\n")
}

#[tauri::command]
fn read_status(run_directory: String) -> Result<Value, String> {
    let path = PathBuf::from(run_directory)
        .join("control")
        .join("status.json");
    let text = fs::read_to_string(path).map_err(|e| e.to_string())?;
    serde_json::from_str(&text).map_err(|e| e.to_string())
}

#[tauri::command]
fn read_run_manifest(run_directory: String) -> Result<Value, String> {
    let path = PathBuf::from(run_directory).join("run.json");
    let text = fs::read_to_string(path).map_err(|e| e.to_string())?;
    serde_json::from_str(&text).map_err(|e| e.to_string())
}

fn viewer_log_tail(path: &Path) -> String {
    let Ok(text) = fs::read_to_string(path) else {
        return String::new();
    };
    let mut characters: Vec<char> = text.chars().rev().take(4000).collect();
    characters.reverse();
    characters
        .into_iter()
        .collect::<String>()
        .trim()
        .to_string()
}

fn wait_for_viewer_server(
    child: &mut Child,
    port: u16,
    log_path: &Path,
    timeout: Duration,
) -> Result<(), String> {
    let address = SocketAddr::V4(SocketAddrV4::new(Ipv4Addr::LOCALHOST, port));
    let deadline = Instant::now() + timeout;
    loop {
        if let Some(status) = child.try_wait().map_err(|e| e.to_string())? {
            let details = viewer_log_tail(log_path);
            return Err(if details.is_empty() {
                format!("3D viewer exited before startup ({status})")
            } else {
                format!("3D viewer exited before startup ({status}):\n{details}")
            });
        }
        if TcpStream::connect_timeout(&address, Duration::from_millis(100)).is_ok() {
            return Ok(());
        }
        if Instant::now() >= deadline {
            let _ = child.kill();
            let _ = child.wait();
            let details = viewer_log_tail(log_path);
            return Err(if details.is_empty() {
                "timed out waiting for the 3D viewer to start".to_string()
            } else {
                format!("timed out waiting for the 3D viewer to start:\n{details}")
            });
        }
        std::thread::sleep(Duration::from_millis(50));
    }
}

#[tauri::command]
fn viewer_heartbeat(
    processes: tauri::State<'_, ViewerProcesses>,
    run_directory: String,
) -> Result<bool, String> {
    let mut children = processes.0.lock().map_err(|e| e.to_string())?;
    let running = if let Some(child) = children.get_mut(&run_directory) {
        child.try_wait().map_err(|e| e.to_string())?.is_none()
    } else {
        false
    };
    if !running {
        children.remove(&run_directory);
        let _ = fs::remove_file(
            PathBuf::from(run_directory)
                .join("control")
                .join("viewer.attached"),
        );
        return Ok(false);
    }
    drop(children);
    atomic_write(
        &PathBuf::from(run_directory)
            .join("control")
            .join("viewer.attached"),
        b"1\n",
    )?;
    Ok(true)
}

#[tauri::command]
fn start_viewer(
    processes: tauri::State<'_, ViewerProcesses>,
    pvpython: String,
    viewer_script: String,
    pythonpath: String,
    run_directory: String,
    language: String,
) -> Result<u16, String> {
    if !Path::new(&pvpython).is_file() {
        return Err(format!("pvpython not found: {pvpython}"));
    }
    if !Path::new(&viewer_script).is_file() {
        return Err(format!("viewer script not found: {viewer_script}"));
    }
    if !matches!(language.as_str(), "en" | "zh-CN") {
        return Err(format!("unsupported viewer language: {language}"));
    }
    let listener = TcpListener::bind("127.0.0.1:0").map_err(|e| e.to_string())?;
    let port = listener.local_addr().map_err(|e| e.to_string())?.port();
    drop(listener);

    let key = run_directory.clone();
    let mut children = processes.0.lock().map_err(|e| e.to_string())?;
    for (old_run, old) in children.iter_mut() {
        let _ = old.kill();
        let _ = old.wait();
        let _ = fs::remove_file(
            PathBuf::from(old_run)
                .join("control")
                .join("viewer.attached"),
        );
    }
    children.clear();
    let log_path = PathBuf::from(&run_directory)
        .join("control")
        .join("viewer.log");
    if let Some(parent) = log_path.parent() {
        fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    let log = fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(&log_path)
        .map_err(|e| e.to_string())?;
    let log_error = log.try_clone().map_err(|e| e.to_string())?;
    let mut child = Command::new(pvpython)
        .arg(viewer_script)
        .arg(&run_directory)
        .arg("--language")
        .arg(language)
        .arg("--server")
        .arg("--port")
        .arg(port.to_string())
        .env("PYTHONPATH", pythonpath)
        .stdout(Stdio::from(log))
        .stderr(Stdio::from(log_error))
        .spawn()
        .map_err(|e| e.to_string())?;
    wait_for_viewer_server(&mut child, port, &log_path, Duration::from_secs(60))?;
    children.insert(key, child);
    drop(children);
    atomic_write(
        &PathBuf::from(run_directory)
            .join("control")
            .join("viewer.attached"),
        b"1\n",
    )?;
    Ok(port)
}

pub fn run() {
    tauri::Builder::default()
        .manage(ViewerProcesses::default())
        .invoke_handler(tauri::generate_handler![
            choose_yaml_file,
            choose_run_directory,
            load_yaml,
            save_yaml,
            validate_config,
            start_simulation,
            send_control,
            read_status,
            read_run_manifest,
            viewer_heartbeat,
            start_viewer,
        ])
        .run(tauri::generate_context!())
        .expect("error while running ATCG3D Studio");
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temporary_directory() -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock")
            .as_nanos();
        std::env::temp_dir().join(format!("atcg3d-studio-test-{}-{nonce}", std::process::id()))
    }

    #[test]
    fn yaml_round_trip_and_control_request_are_atomic() {
        let directory = temporary_directory();
        let yaml_path = directory.join("config.yaml");
        let value = serde_json::json!({
            "schema": "atcg3d.model_config",
            "simulation": {"threads": 18},
            "output": {"enabled": true}
        });
        save_yaml(yaml_path.to_string_lossy().into_owned(), value.clone()).expect("save YAML");
        assert_eq!(
            load_yaml(yaml_path.to_string_lossy().into_owned()).expect("load YAML"),
            value
        );
        assert!(!yaml_path.with_extension("yaml.tmp").exists());

        send_control(
            directory.to_string_lossy().into_owned(),
            "pause".to_string(),
        )
        .expect("write pause request");
        assert_eq!(
            fs::read_to_string(directory.join("control").join("pause.request"))
                .expect("read pause request"),
            "1\n"
        );
        assert!(send_control(
            directory.to_string_lossy().into_owned(),
            "unknown".to_string()
        )
        .is_err());
        fs::remove_dir_all(directory).expect("remove test directory");
    }

    #[test]
    fn viewer_readiness_reports_early_exit() {
        let mut exited = Command::new("/usr/bin/true")
            .spawn()
            .expect("spawn exiting process");
        let error = wait_for_viewer_server(
            &mut exited,
            9,
            Path::new("/nonexistent/viewer.log"),
            Duration::from_secs(1),
        )
        .expect_err("exited viewer should fail readiness");
        assert!(error.contains("exited before startup"));
    }
}

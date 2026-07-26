const invoke = window.__TAURI__.core.invoke;

const translations = {
  en: {
    "language.label": "Language",
    "run.notSelected": "No run directory selected",
    "run.notRunning": "Not running",
    "viewport.title": "3D simulation and visualization",
    "viewport.hint": "Load a YAML configuration or open an existing run from the panel on the right.",
    "viewport.loading": "Starting ParaView and loading the latest 3D frame…",
    "viewport.loadFailed": "The 3D view did not load. Check the run log for details, then retry Open 3D view.",
    "viewport.frameTitle": "ATCG3D 3D viewer",
    "tabs.run": "Run",
    "tabs.parameters": "Parameters",
    "tabs.effective": "Effective",
    "tabs.diagnostics": "Diagnostics",
    "run.filesAndPrograms": "Files and programs",
    "run.configYaml": "Configuration YAML",
    "common.choose": "Choose",
    "run.runDirectory": "Run directory",
    "common.open": "Open",
    "run.programPaths": "Program paths",
    "run.launcherScript": "Detached launcher script",
    "run.viewerScript": "Viewer script",
    "run.loadParameters": "Load parameters",
    "run.validateParameters": "Validate parameters",
    "run.saveYaml": "Save YAML",
    "run.startSimulation": "Start simulation",
    "run.pauseSimulation": "Pause simulation",
    "run.resumeSimulation": "Resume simulation",
    "run.checkpointStop": "Checkpoint and stop",
    "run.terminate": "Stop (no new checkpoint)",
    "run.openViewer": "Open 3D view",
    "run.log": "Run log",
    "parameters.edit": "Edit parameters",
    "parameters.filter": "Filter parameters…",
    "parameters.hint": "All simulation parameters remain in YAML. Objects and arrays can be expanded, and values retain their original types.",
    "effective.title": "Effective configuration",
    "effective.hint": "The configuration strictly validated by C++ is shown here after validation or startup.",
    "effective.notValidated": "Not validated yet.",
    "diagnostics.title": "Run diagnostics",
    "diagnostics.noStatus": "No status available.",
    "message.waiting": "Waiting for an action.",
    "message.arrayInvalid": "Invalid array syntax at {path}.",
    "message.parametersLoaded": "Parameters loaded.",
    "message.loadFirst": "Load the parameters first.",
    "message.yamlSaved": "YAML saved atomically.",
    "message.validationPassed": "Strict C++ validation passed.",
    "message.simulationStarted": "Simulation launcher started.\n{details}",
    "message.controlSubmitted": "Control request submitted: {action}.",
    "message.viewerStarted": "3D viewer started on local port {port}.",
    "message.runOpened": "Existing run opened and its effective configuration loaded.",
    "message.runManifestFailed": "Failed to read the run configuration: {error}",
    "message.viewerRestarting": "Applying the new language to the 3D viewer…",
    "action.pause": "pause",
    "action.resume": "resume",
    "action.checkpoint_stop": "checkpoint and stop",
    "action.terminate": "stop without a new checkpoint",
    "status.state": "state",
    "status.time": "time",
    "status.cells": "cells",
    "status.vessels": "vessels",
    "status.speed": "speed",
    "status.simHoursPerWallHour": "sim h/wall h",
    "status.idle": "idle",
    "status.running": "running",
    "status.paused": "paused",
    "status.stopping": "stopping",
    "status.completed": "completed",
    "status.failed": "failed",
  },
  "zh-CN": {
    "language.label": "语言",
    "run.notSelected": "未选择运行目录",
    "run.notRunning": "未运行",
    "viewport.title": "三维模拟与可视化",
    "viewport.hint": "在右侧加载 YAML 配置或打开已有运行结果。",
    "viewport.loading": "正在启动 ParaView 并加载最新三维帧…",
    "viewport.loadFailed": "三维视图未能加载。请检查运行日志，然后重试“打开三维视图”。",
    "viewport.frameTitle": "ATCG3D 三维查看器",
    "tabs.run": "运行",
    "tabs.parameters": "参数",
    "tabs.effective": "当前参数",
    "tabs.diagnostics": "诊断",
    "run.filesAndPrograms": "文件与程序",
    "run.configYaml": "配置 YAML",
    "common.choose": "选择",
    "run.runDirectory": "运行目录",
    "common.open": "打开",
    "run.programPaths": "程序路径",
    "run.launcherScript": "后台启动脚本",
    "run.viewerScript": "Viewer 脚本",
    "run.loadParameters": "加载参数",
    "run.validateParameters": "验证参数",
    "run.saveYaml": "保存 YAML",
    "run.startSimulation": "开始模拟",
    "run.pauseSimulation": "暂停模拟",
    "run.resumeSimulation": "继续模拟",
    "run.checkpointStop": "Checkpoint 并停止",
    "run.terminate": "停止（无新 checkpoint）",
    "run.openViewer": "打开三维视图",
    "run.log": "运行日志",
    "parameters.edit": "参数编辑",
    "parameters.filter": "筛选参数…",
    "parameters.hint": "所有模拟参数仍写入 YAML。对象和数组可展开，数值保持原类型。",
    "effective.title": "当前生效配置",
    "effective.hint": "验证或开始运行后，这里显示由 C++ 严格验证的实际参数。",
    "effective.notValidated": "尚未验证。",
    "diagnostics.title": "运行诊断",
    "diagnostics.noStatus": "尚无状态。",
    "message.waiting": "等待操作。",
    "message.arrayInvalid": "数组格式错误：{path}。",
    "message.parametersLoaded": "参数已加载。",
    "message.loadFirst": "请先加载参数。",
    "message.yamlSaved": "YAML 已原子保存。",
    "message.validationPassed": "C++ 严格验证通过。",
    "message.simulationStarted": "模拟启动器已启动。\n{details}",
    "message.controlSubmitted": "已提交控制命令：{action}。",
    "message.viewerStarted": "三维 Viewer 已在本地端口 {port} 启动。",
    "message.runOpened": "已打开现有运行目录并读取当前参数。",
    "message.runManifestFailed": "运行目录参数读取失败：{error}",
    "message.viewerRestarting": "正在将新语言应用到三维 Viewer…",
    "action.pause": "暂停",
    "action.resume": "继续",
    "action.checkpoint_stop": "Checkpoint 并停止",
    "action.terminate": "停止且不创建新 checkpoint",
    "status.state": "状态",
    "status.time": "时间",
    "status.cells": "细胞",
    "status.vessels": "血管",
    "status.speed": "速度",
    "status.simHoursPerWallHour": "模拟小时/现实小时",
    "status.idle": "空闲",
    "status.running": "运行中",
    "status.paused": "已暂停",
    "status.stopping": "正在停止",
    "status.completed": "已完成",
    "status.failed": "失败",
  },
};

const supportedLanguages = new Set(Object.keys(translations));
let currentLanguage = localStorage.getItem("language");
if (!supportedLanguages.has(currentLanguage)) currentLanguage = "en";

function t(key, parameters = {}) {
  const template = translations[currentLanguage][key]
    ?? translations.en[key] ?? key;
  return template.replace(/\{(\w+)\}/g, (_, name) =>
    Object.prototype.hasOwnProperty.call(parameters, name)
      ? String(parameters[name]) : `{${name}}`);
}

const defaults = {
  binary: "/Volumes/Work_Active/simulation/ver7/build/atcg3d",
  launcher: "/Users/taolee/Documents/GitHub/Ver7/scripts/launch_atcg3d_detached.sh",
  config: "/Users/taolee/Documents/GitHub/Ver7/configs/single_r_stage0_2160h_seed1.yaml",
  run: "/Volumes/Work_Active/simulation/ver7/run_single_r_stage0_2160h_seed1",
  pvpython: "/Applications/ParaView-6.1.1.app/Contents/bin/pvpython",
  viewer: "/Users/taolee/Documents/GitHub/Ver7/visualization/viewer/app.py",
  pythonpath: `${localStorage.getItem("HOME") || "/Users/taolee"}/.pyenv/versions/atcg3d-paraview-3.12.7/lib/python3.12/site-packages`,
};

let configObject = null;
let effectiveObject = null;
let viewerActive = false;
let viewerLoading = false;
let viewerLoadFailed = false;
let latestStatus = null;
let lastMessage = { key: "message.waiting", parameters: {}, raw: null, isError: false };

const $ = (id) => document.getElementById(id);
const fields = {
  binary: $("binary-path"), launcher: $("launcher-path"),
  config: $("config-path"), run: $("run-directory"),
  pvpython: $("pvpython-path"), viewer: $("viewer-path"),
  pythonpath: $("pythonpath"),
};
for (const [key, input] of Object.entries(fields)) {
  input.value = localStorage.getItem(key) || defaults[key];
  input.addEventListener("change", () => localStorage.setItem(key, input.value));
}

function renderMessage() {
  $("message").textContent = lastMessage.raw === null
    ? t(lastMessage.key, lastMessage.parameters)
    : String(lastMessage.raw);
  $("message").style.color = lastMessage.isError ? "#ff9b9b" : "";
}

function messageKey(key, parameters = {}, isError = false) {
  lastMessage = { key, parameters, raw: null, isError };
  renderMessage();
}

function messageError(error, prefixKey = null) {
  const detail = String(error);
  lastMessage = prefixKey
    ? { key: prefixKey, parameters: { error: detail }, raw: null, isError: true }
    : { key: "", parameters: {}, raw: detail, isError: true };
  renderMessage();
}

function localizedState(value) {
  const normalized = String(value || "idle").toLowerCase();
  const key = `status.${normalized}`;
  return translations[currentLanguage][key]
    ?? translations.en[key] ?? String(value || "idle");
}

function renderStatus(status = latestStatus) {
  if (!status) {
    $("state").textContent = `${t("status.state")}: ${t("status.idle")}`;
    $("time").textContent = `${t("status.time")}=0 h`;
    $("cells").textContent = `${t("status.cells")}=0`;
    $("vessels").textContent = `${t("status.vessels")}=0`;
    $("speed").textContent =
      `${t("status.speed")}=0 ${t("status.simHoursPerWallHour")}`;
    $("eta").textContent = "ETA=—";
    $("progress-label").textContent = t("run.notRunning");
    $("run-name").textContent = t("run.notSelected");
    return;
  }
  $("state").textContent =
    `${t("status.state")}: ${localizedState(status.state)}`;
  $("time").textContent =
    `${t("status.time")}=${Number(status.time_hours).toFixed(3)} h`;
  $("cells").textContent =
    `${t("status.cells")}=${Number(status.alive_cells).toLocaleString(currentLanguage)}`;
  $("vessels").textContent =
    `${t("status.vessels")}=${Number(status.vessel_nodes).toLocaleString(currentLanguage)}`;
  $("speed").textContent =
    `${t("status.speed")}=${Number(status.simulated_hours_per_wall_hour).toFixed(2)} ${t("status.simHoursPerWallHour")}`;
  $("eta").textContent = status.eta_wall_seconds >= 0
    ? `ETA=${(status.eta_wall_seconds / 3600).toFixed(1)} h` : "ETA=—";
  $("progress").value = status.progress_fraction || 0;
  $("progress-label").textContent =
    `${(100 * (status.progress_fraction || 0)).toFixed(2)}%`;
  $("run-name").textContent =
    fields.run.value.split("/").filter(Boolean).pop() || t("run.notSelected");
  $("diagnostics").textContent = JSON.stringify(status, null, 2);
}

function renderEmptyView() {
  $("empty-view-title").textContent = t("viewport.title");
  $("empty-view-hint").textContent = viewerLoading
    ? t("viewport.loading")
    : viewerLoadFailed ? t("viewport.loadFailed") : t("viewport.hint");
}

function applyLanguage() {
  document.documentElement.lang = currentLanguage;
  $("language-select").value = currentLanguage;
  $("language-select").setAttribute("aria-label", t("language.label"));
  document.querySelectorAll("[data-i18n]").forEach(element => {
    element.textContent = t(element.dataset.i18n);
  });
  document.querySelectorAll("[data-i18n-placeholder]").forEach(element => {
    element.placeholder = t(element.dataset.i18nPlaceholder);
  });
  document.querySelectorAll("[data-i18n-title]").forEach(element => {
    element.title = t(element.dataset.i18nTitle);
  });
  if (effectiveObject) {
    $("effective-config").textContent = JSON.stringify(effectiveObject, null, 2);
  }
  renderEmptyView();
  renderStatus();
  renderMessage();
}

function updatePathFromConfig() {
  const run = configObject?.output?.directory;
  if (typeof run === "string" && run) {
    fields.run.value = run;
    localStorage.setItem("run", run);
  }
}

function createScalarEditor(parent, key, value, path) {
  const row = document.createElement("div");
  row.className = "tree-row";
  row.dataset.path = path.toLowerCase();
  const label = document.createElement("label");
  label.textContent = key;
  const input = document.createElement("input");
  if (typeof value === "boolean") {
    input.type = "checkbox";
    input.checked = value;
    input.addEventListener("change", () => parent[key] = input.checked);
  } else {
    input.type = typeof value === "number" ? "number" : "text";
    if (typeof value === "number") input.step = "any";
    input.value = value === null ? "" : String(value);
    input.addEventListener("change", () => {
      if (typeof value === "number") {
        const parsed = Number(input.value);
        parent[key] = Number.isFinite(parsed) ? parsed : input.value;
      } else if (value === null && input.value === "") {
        parent[key] = null;
      } else {
        parent[key] = input.value;
      }
      updatePathFromConfig();
    });
  }
  row.append(label, input);
  return row;
}

function renderObject(container, object, path = "") {
  for (const [key, value] of Object.entries(object || {})) {
    const current = path ? `${path}.${key}` : key;
    if (value && typeof value === "object" && !Array.isArray(value)) {
      const details = document.createElement("details");
      details.open = path === "";
      details.dataset.path = current.toLowerCase();
      const summary = document.createElement("summary");
      summary.className = "tree-key";
      summary.textContent = key;
      const group = document.createElement("div");
      group.className = "tree-group";
      renderObject(group, value, current);
      details.append(summary, group);
      container.append(details);
    } else if (Array.isArray(value)) {
      const row = createScalarEditor(object, key, JSON.stringify(value), current);
      const input = row.querySelector("input");
      input.addEventListener("change", () => {
        try { object[key] = JSON.parse(input.value); }
        catch { messageKey("message.arrayInvalid", { path: current }, true); }
      });
      container.append(row);
    } else {
      container.append(createScalarEditor(object, key, value, current));
    }
  }
}

function renderParameters() {
  const tree = $("parameter-tree");
  tree.innerHTML = "";
  renderObject(tree, configObject);
}

async function loadConfig() {
  try {
    configObject = await invoke("load_yaml", { path: fields.config.value });
    renderParameters();
    updatePathFromConfig();
    messageKey("message.parametersLoaded");
  } catch (error) { messageError(error); }
}

async function saveConfig() {
  if (!configObject) {
    messageKey("message.loadFirst", {}, true);
    return false;
  }
  try {
    if (configObject.output && typeof configObject.output === "object") {
      configObject.output.directory = fields.run.value;
    }
    await invoke("save_yaml", { path: fields.config.value, value: configObject });
    messageKey("message.yamlSaved");
    return true;
  } catch (error) {
    messageError(error);
    return false;
  }
}

async function validateConfig() {
  if (!(await saveConfig())) return false;
  try {
    const text = await invoke("validate_config", {
      binary: fields.binary.value, config: fields.config.value,
    });
    effectiveObject = JSON.parse(text);
    $("effective-config").textContent = JSON.stringify(effectiveObject, null, 2);
    messageKey("message.validationPassed");
    return true;
  } catch (error) {
    messageError(error);
    return false;
  }
}

async function startSimulation() {
  if (!(await validateConfig())) return;
  try {
    const result = await invoke("start_simulation", {
      launcher: fields.launcher.value, binary: fields.binary.value,
      config: fields.config.value, runDirectory: fields.run.value,
    });
    messageKey("message.simulationStarted", { details: result });
    setTimeout(openViewer, 2500);
  } catch (error) { messageError(error); }
}

async function control(action) {
  try {
    await invoke("send_control", {
      runDirectory: fields.run.value, action,
    });
    messageKey("message.controlSubmitted", {
      action: t(`action.${action}`),
    });
  } catch (error) { messageError(error); }
}

async function openViewer() {
  const frame = $("viewer");
  viewerActive = false;
  viewerLoading = true;
  viewerLoadFailed = false;
  frame.onload = null;
  frame.onerror = null;
  frame.style.display = "none";
  frame.src = "about:blank";
  $("empty-view").style.display = "grid";
  renderEmptyView();
  try {
    const port = await invoke("start_viewer", {
      pvpython: fields.pvpython.value,
      viewerScript: fields.viewer.value,
      pythonpath: fields.pythonpath.value,
      runDirectory: fields.run.value,
      language: currentLanguage,
    });
    viewerActive = true;
    let completed = false;
    const finishLoading = () => {
      if (completed) return;
      completed = true;
      viewerLoading = false;
      viewerLoadFailed = false;
      frame.style.display = "block";
      $("empty-view").style.display = "none";
      messageKey("message.viewerStarted", { port });
    };
    frame.onload = finishLoading;
    frame.onerror = () => {
      if (completed) return;
      completed = true;
      viewerLoading = false;
      viewerLoadFailed = true;
      renderEmptyView();
    };
    frame.src = `http://127.0.0.1:${port}/?studio=${Date.now()}`;
    setTimeout(() => {
      if (completed) return;
      completed = true;
      viewerLoading = false;
      viewerLoadFailed = true;
      renderEmptyView();
    }, 15000);
  } catch (error) {
    viewerLoading = false;
    viewerLoadFailed = true;
    renderEmptyView();
    messageError(error);
  }
}

async function pollStatus() {
  if (viewerActive) {
    invoke("viewer_heartbeat", {
      runDirectory: fields.run.value,
    }).then(running => { viewerActive = Boolean(running); })
      .catch(() => { viewerActive = false; });
  }
  try {
    const status = await invoke("read_status", { runDirectory: fields.run.value });
    latestStatus = status;
    renderStatus(status);
  } catch (_) {}
}

document.querySelectorAll(".tabs button").forEach(button => {
  button.addEventListener("click", () => {
    document.querySelectorAll(".tabs button").forEach(v => v.classList.remove("active"));
    document.querySelectorAll(".tab").forEach(v => v.classList.remove("active"));
    button.classList.add("active");
    $(`tab-${button.dataset.tab}`).classList.add("active");
  });
});
$("parameter-filter").addEventListener("input", event => {
  const query = event.target.value.trim().toLowerCase();
  document.querySelectorAll("#parameter-tree [data-path]").forEach(node => {
    node.classList.toggle("hidden", query && !node.dataset.path.includes(query));
  });
});
$("choose-config").onclick = async () => {
  const path = await invoke("choose_yaml_file");
  if (path) { fields.config.value = path; await loadConfig(); }
};
$("choose-run").onclick = async () => {
  const path = await invoke("choose_run_directory");
  if (path) {
    fields.run.value = path;
    localStorage.setItem("run", path);
    try {
      const manifest = await invoke("read_run_manifest", {
        runDirectory: path,
      });
      effectiveObject = manifest.effective_config || manifest;
      $("effective-config").textContent =
        JSON.stringify(effectiveObject, null, 2);
      messageKey("message.runOpened");
    } catch (error) {
      messageError(error, "message.runManifestFailed");
    }
    await openViewer();
  }
};
$("load-config").onclick = loadConfig;
$("save-config").onclick = saveConfig;
$("validate").onclick = validateConfig;
$("start").onclick = startSimulation;
$("pause").onclick = () => control("pause");
$("resume").onclick = () => control("resume");
$("checkpoint-stop").onclick = () => control("checkpoint_stop");
$("terminate").onclick = () => control("terminate");
$("open-viewer").onclick = openViewer;
$("language-select").addEventListener("change", async event => {
  const requested = event.target.value;
  if (!supportedLanguages.has(requested) || requested === currentLanguage) return;
  currentLanguage = requested;
  localStorage.setItem("language", currentLanguage);
  applyLanguage();
  if (viewerActive) {
    messageKey("message.viewerRestarting");
    await openViewer();
  }
});

async function initialize() {
  applyLanguage();
  await loadConfig();
  await pollStatus();
  if (latestStatus && !viewerActive && !viewerLoading) {
    await openViewer();
  }
}

initialize();
setInterval(pollStatus, 1000);

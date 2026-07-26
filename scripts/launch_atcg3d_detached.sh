#!/bin/sh
set -eu

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 ATCG3D_BINARY CONFIG_YAML CONTROL_PREFIX" >&2
    echo "Creates CONTROL_PREFIX.pid, .log, and atomic .exit files." >&2
    exit 2
fi

binary=$1
config=$2
prefix=$3
pid_file="${prefix}.pid"
log_file="${prefix}.log"
exit_file="${prefix}.exit"

if [ ! -x "$binary" ]; then
    echo "ATCG3D binary is not executable: $binary" >&2
    exit 2
fi
if [ ! -f "$config" ]; then
    echo "Configuration file does not exist: $config" >&2
    exit 2
fi
if [ -f "$pid_file" ]; then
    old_pid=$(sed -n '1p' "$pid_file")
    if [ -n "$old_pid" ] && kill -0 "$old_pid" 2>/dev/null; then
        echo "A managed process is already running with PID $old_pid" >&2
        exit 3
    fi
fi

mkdir -p "$(dirname "$prefix")"
rm -f "$exit_file" "${exit_file}.tmp" "${pid_file}.tmp"

nohup /bin/sh -c '
    binary=$1
    config=$2
    exit_file=$3
    "$binary" --config "$config"
    code=$?
    {
        printf "exit_code=%s\n" "$code"
        printf "finished_utc=%s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    } > "${exit_file}.tmp"
    mv "${exit_file}.tmp" "$exit_file"
    exit "$code"
' atcg3d-launcher "$binary" "$config" "$exit_file" >> "$log_file" 2>&1 </dev/null &

launcher_pid=$!
printf "%s\n" "$launcher_pid" > "${pid_file}.tmp"
mv "${pid_file}.tmp" "$pid_file"
echo "ATCG3D launcher PID: $launcher_pid"
echo "Log: $log_file"
echo "Exit status: $exit_file"

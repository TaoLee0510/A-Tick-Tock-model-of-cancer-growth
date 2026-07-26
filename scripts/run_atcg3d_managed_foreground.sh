#!/bin/sh
set -u

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 ATCG3D_BINARY CONFIG_YAML CONTROL_PREFIX" >&2
    exit 2
fi

binary=$1
config=$2
prefix=$3
pid_file="${prefix}.pid"
exit_file="${prefix}.exit"
child=
finished=0

write_exit() {
    code=$1
    if [ "$finished" -eq 1 ]; then
        return
    fi
    finished=1
    {
        printf "exit_code=%s\n" "$code"
        printf "finished_utc=%s\n" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    } > "${exit_file}.tmp"
    mv "${exit_file}.tmp" "$exit_file"
}

terminate_child() {
    if [ -n "$child" ] && kill -0 "$child" 2>/dev/null; then
        kill -TERM "$child" 2>/dev/null || true
        wait "$child" 2>/dev/null || true
    fi
    write_exit 143
    exit 143
}

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

"$binary" --config "$config" &
child=$!
printf "%s\n" "$child" > "${pid_file}.tmp"
mv "${pid_file}.tmp" "$pid_file"

trap terminate_child HUP INT TERM

wait "$child"
code=$?
write_exit "$code"
exit "$code"

#!/bin/sh

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 ATCG3D_BINARY CONFIG_YAML CONTROL_PREFIX" >&2
    exit 2
fi

binary=$1
config=$2
prefix=$3
pid_file="${prefix}.pid"
exit_file="${prefix}.exit"
child_pid=
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

terminate() {
    if [ -n "$child_pid" ] && kill -0 "$child_pid" 2>/dev/null; then
        kill -TERM "$child_pid" 2>/dev/null || true
        wait "$child_pid" 2>/dev/null || true
    fi
    write_exit 143
    exit 143
}

trap terminate HUP INT TERM
rm -f "$exit_file" "${exit_file}.tmp" "${pid_file}.tmp"

"$binary" --config "$config" &
child_pid=$!
printf "%s\n" "$child_pid" > "${pid_file}.tmp"
mv "${pid_file}.tmp" "$pid_file"

wait "$child_pid"
code=$?
write_exit "$code"
exit "$code"

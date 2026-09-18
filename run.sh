#!/bin/sh
# Cross-platform (Linux/macOS) equivalent of run.bat.
echo "=== Running LYNXgwas ==="

CP="bin"
for f in lib/*.jar; do
  [ -e "$f" ] && CP="$CP:$f"
done

# run.bat opens the browser only after the java process exits, since a plain blocking
# `java ... Main` call never returns control to the batch script while the server is up —
# it likely never fires at the intended moment there. Backgrounding it here lets us open the
# browser once the server is actually likely to be listening, matching what the PyPI launcher
# (python/lynxgwas/cli.py) already does correctly.
java -Xmx8g -cp "$CP" Main "$@" &
JAVA_PID=$!
sleep 2
if kill -0 "$JAVA_PID" 2>/dev/null; then
  echo
  echo "Opening home page..."
  if command -v xdg-open >/dev/null 2>&1; then
    xdg-open http://localhost:8765/ >/dev/null 2>&1
  elif command -v open >/dev/null 2>&1; then
    open http://localhost:8765/
  else
    echo "Open http://localhost:8765/ in your browser."
  fi
fi

wait "$JAVA_PID"
status=$?
if [ $status -ne 0 ]; then
  echo "[FAIL] Run failed. See above for errors."
  exit 1
fi

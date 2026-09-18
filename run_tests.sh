#!/bin/sh
# Cross-platform (Linux/macOS) equivalent of run_tests.bat.
set -e
echo "=== Running LYNXgwas tests ==="

LIBCP=""
for f in lib/*.jar; do
  [ -e "$f" ] && LIBCP="$LIBCP:$f"
done

javac -encoding UTF-8 -d bin -cp "bin$LIBCP" \
  tests/MultiLocusScannerTest.java \
  tests/LdCalculatorTest.java \
  tests/GwasQcTest.java

FAILED=0

echo
echo "--- MultiLocusScannerTest ---"
java -cp "bin$LIBCP" MultiLocusScannerTest || FAILED=1

echo
echo "--- LdCalculatorTest ---"
java -cp "bin$LIBCP" LdCalculatorTest || FAILED=1

echo
echo "--- GwasQcTest ---"
java -cp "bin$LIBCP" GwasQcTest || FAILED=1

echo
if [ "$FAILED" -eq 0 ]; then
  echo "[OK] All test suites passed."
else
  echo "[FAIL] One or more test suites failed."
  exit 1
fi

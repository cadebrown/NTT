#!/bin/sh


# how many hex digits
#HEXDIGS=$((32))
if [ -z "${HEXDIGS}" ]; then
    HEXDIGS=$((1024))
fi

rand() {
    openssl rand -hex $HEXDIGS
}

A=$(rand)
B=$(rand)

echo $A > /tmp/A.txt
echo $B > /tmp/B.txt

./tools/mul_py.py /tmp/A.txt /tmp/B.txt > /tmp/C_py.txt
./bin/ntt mulhex /tmp/A.txt /tmp/B.txt > /tmp/C_ntt.txt

# ensure they are the same output
cmp /tmp/C_py.txt /tmp/C_ntt.txt && echo "Success!" || echo "Failure!"
echo "Check /tmp/A.txt /tmp/B.txt /tmp/C_py.txt /tmp/C_ntt.txt"
echo "Run with 'HEXDIGS=1234 $0' to test different sizes"

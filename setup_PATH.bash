# source this file to add Mhp1 scripts to PATH
BINDIR=${HOME}/Projects/Beckstein/Mhp1/Library/bin

test -d ${BINDIR} || { echo "BINDIR=${BINDIR} not found."; return 1; }

# better: only add to PATH if not there yet ...
# but I am too lazy right now

export PATH=${BINDIR}:${PATH}

echo "-- prepended ${BINDIR} to PATH"


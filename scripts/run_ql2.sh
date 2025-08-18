#!/bin/bash
# Source the user's shell configuration to get aliases and PATH
source ~/.bashrc 2>/dev/null || source ~/.bash_profile 2>/dev/null || true

# Try to find IDL in common locations or use PATH
if command -v idl >/dev/null 2>&1; then
     IDL_CMD="idl"
# elif [ -f "/Applications/harris/idl87/bin/idl" ]; then
#     IDL_CMD="/Applications/harris/idl87/bin/idl"
# elif [ -f "/Applications/harris/idl88/bin/idl" ]; then
#     IDL_CMD="/Applications/harris/idl88/bin/idl"
# elif [ -f "/Applications/harris/idl89/bin/idl" ]; then
#     IDL_CMD="/Applications/harris/idl89/bin/idl"
else
     echo "Error: IDL not found. Please ensure IDL is installed and in your PATH."
     exit 1
fi

${IDL_CMD} ${OSIRIS_ROOT}/ql2/qlook2_startup.pro

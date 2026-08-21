#!/bin/bash
set -e

cleanup() {
    echo "Received shutdown signal, cleaning up..."
    if [ ! -z "$SERVER_PID" ]; then
        kill -TERM "$SERVER_PID" 2>/dev/null || true
        wait "$SERVER_PID" 2>/dev/null || true
    fi
    exit 0
}
trap cleanup SIGTERM SIGINT

MCP_HOST=${MCP_HOST:-"0.0.0.0"}
MCP_PORT=${MCP_PORT:-"8080"}

export PYTHONPATH="/app:${PYTHONPATH}"

echo "Checking dependencies..."
python -c "import mcp; print('mcp SDK available')"
python -c "from rdkit import Chem; print('RDKit available')"
echo "All dependencies available"

echo "Starting RDKit Tools MCP Server..."
echo "Host: $MCP_HOST"
echo "Port: $MCP_PORT"

python server.py --host "$MCP_HOST" --port "$MCP_PORT" &
SERVER_PID=$!

sleep 3

if kill -0 "$SERVER_PID" 2>/dev/null; then
    echo "Server started successfully (PID: $SERVER_PID)"
    wait "$SERVER_PID"
else
    echo "Failed to start server"
    exit 1
fi

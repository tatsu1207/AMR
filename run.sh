#!/bin/bash
cd "$(dirname "$0")"
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate amr

case "${1:-start}" in
    start)
        if pgrep -f "uvicorn app.main:app.*--port 8000" > /dev/null; then
            echo "AMR Predictor is already running (PID $(pgrep -f 'uvicorn app.main:app.*--port 8000'))"
            exit 1
        fi
        nohup uvicorn app.main:app --host 0.0.0.0 --port 8000 > uvicorn.log 2>&1 &
        echo "AMR Predictor started (PID $!) on http://0.0.0.0:8000"
        echo "Logs: uvicorn.log"
        ;;
    stop)
        pkill -f "uvicorn app.main:app.*--port 8000" && echo "AMR Predictor stopped" || echo "Not running"
        ;;
    restart)
        pkill -f "uvicorn app.main:app.*--port 8000" 2>/dev/null && echo "Stopped old process" && sleep 1
        nohup uvicorn app.main:app --host 0.0.0.0 --port 8000 > uvicorn.log 2>&1 &
        echo "AMR Predictor restarted (PID $!) on http://0.0.0.0:8000"
        echo "Logs: uvicorn.log"
        ;;
    status)
        if pgrep -f "uvicorn app.main:app.*--port 8000" > /dev/null; then
            echo "AMR Predictor is running (PID $(pgrep -f 'uvicorn app.main:app.*--port 8000'))"
        else
            echo "AMR Predictor is not running"
        fi
        ;;
    *)
        echo "Usage: $0 {start|stop|restart|status}"
        exit 1
        ;;
esac

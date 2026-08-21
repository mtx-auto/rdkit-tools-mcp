FROM python:3.11-slim

WORKDIR /app

ENV DEBIAN_FRONTEND=noninteractive

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY server.py .
COPY entrypoint.sh .
COPY README.md .

RUN chmod +x entrypoint.sh

EXPOSE 8080

ENV MCP_HOST=0.0.0.0
ENV MCP_PORT=8080

ENTRYPOINT ["./entrypoint.sh"]

FROM python:3.13-slim

RUN apt-get update && apt-get install -y \
    samtools \
    && rm -rf /var/lib/apt/lists/*

ADD . /antigenfinder
RUN pip install --no-cache-dir ./antigenfinder

ENTRYPOINT ["python3", "-m", "antigenfinder"]
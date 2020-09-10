FROM gcr.io/broad-getzlab-workflows/base_image:latest

WORKDIR build
RUN git clone https://github.com/getzlab/CApy.git && pip install -e ./CApy

WORKDIR /app
ENV PATH=$PATH:/app
COPY hetpull.py .

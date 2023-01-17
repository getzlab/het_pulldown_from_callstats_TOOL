FROM gcr.io/broad-getzlab-workflows/base_image:latest

WORKDIR build
RUN git clone https://github.com/getzlab/CApy.git && cd CApy && git checkout eb95ce4 && pip install .

WORKDIR /app
ENV PATH=$PATH:/app
COPY hetpull.py .

FROM python:3.6.3

RUN apt-get update -y
RUN apt-get install -y python-pip python-dev build-essential
RUN pip install --upgrade pip
RUN apt-get install -y nano

WORKDIR notebooks

COPY python/exampleAnalysis.ipynb ./notebooks/exampleAnalysis.ipynb
COPY python/samiPostAnalysis.py ./notebooks/samiPostAnalysis.py
COPY analysis ./analysis

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# ["mkdir", "notebooks"]
#COPY conf/.jupyter /root/.jupyter
#COPY run_jupyter.sh /

# Jupyter and Tensorboard ports
#EXPOSE 8888 6006

# Store notebooks in this mounted directory
VOLUME notebooks

#CMD ["/run_jupyter.sh"]

CMD ["jupyter", "notebook", "--port=8080", "--no-browser", "--ip=0.0.0.0", "--allow-root"]

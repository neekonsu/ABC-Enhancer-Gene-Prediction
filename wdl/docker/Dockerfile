FROM python:3.7

WORKDIR /usr/src/app

# Setup the python requirements
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# Copy the required scripts
COPY src/ src/

# Setup samtools
RUN wget -O samtools-0.1.19.tar.bz2 https://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2/download
RUN tar xjf samtools-0.1.19.tar.bz2
WORKDIR /usr/src/app/samtools-0.1.19.tar.bz2
RUN make
ENV PATH=/usr/src/app/samtools-0.1.19/:${PATH}
WORKDIR /usr/src/app/


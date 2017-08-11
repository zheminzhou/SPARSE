# apt packages
sudo apt-get install -y make build-essential libssl-dev zlib1g-dev libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncurses5-dev libncursesw5-dev xz-utils tk-dev 

# pyenv & pyenv-virtualenv
git clone https://github.com/yyuu/pyenv.git ~/.pyenv && git clone https://github.com/yyuu/pyenv-virtualenv.git ~/.pyenv/plugins/pyenv-virtualenv
(echo 'export PYENV_ROOT="$HOME/.pyenv"'; echo 'export PATH="$PYENV_ROOT/bin:$PATH"'; echo 'eval "$(pyenv init -)"'; echo 'eval "$(pyenv virtualenv-init -)"' ) >> ~/.bashrc
exec "$SHELL"
pyenv install 2.7.9 && pyenv virtualenv 2.7.9 aDNA
pyenv activate aDNA && pyenv rehash

# pip
pip install cython numpy scipy pandas ujson msgpack-python

# dependencies
## samtools
sudo apt-get install liblzma-dev 
wget https://downloads.sourceforge.net/project/samtools/samtools/1.5/samtools-1.5.tar.bz2 && tar xfj samtools-1.5.tar.bz2  && (cd samtools-1.5 && make)

## bowtie2
sudo apt-get install libtbb-dev
wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.2/bowtie2-2.3.2-source.zip && unzip bowtie2-2.3.2-source.zip && (cd bowtie2-2.3.2 && make)

## capnproto
sudo apt-get install autoconf automake libtool
git clone https://github.com/sandstorm-io/capnproto.git && (cd capnproto/c++ && autoreconf -i && ./configure && make -j6 check && sudo make install)
git clone https://github.com/jparyani/pycapnp.git && pip install --install-option '--force-cython' ./pycapnp

## mash
git clone https://github.com/marbl/Mash.git && (cd Mash && ./bootstrap.sh && ./configure && make && sudo make install)

## SPARSE
git clone https://github.com/zheminzhou/SPARSE.git
pyenv activate aDNA
python 01_db_create.py dbname=refseq
python 02_db_fill.py dbname=refseq update=True

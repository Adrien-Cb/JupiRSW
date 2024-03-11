# JupiRSW

This is a rotating shallow water model for studying vortices at a planet's pole.


https://github.com/Adrien-Cb/JupiRSW/assets/130017358/4483151a-a1b7-42d5-8c36-23c77cf84599


*Work in progress*

## Install

### Libraries 

JupiRSW uses **NumPy**, **Xarray**, **SciPy** and **Matplotlib**.

### ffmpeg

**ffmpeg** is required to save animations as MP4. 
There are many builds available, you can for instance get one [here](https://www.gyan.dev/ffmpeg/builds/).
Once installed, add the path to the ffmpeg executable in `config.py`:
```python
FFMPEG_PATH = R'C:\ffmpeg\bin\ffmpeg.exe'
```

## References and credits

The RK3 scheme implementation is taken from [PyRSW](https://github.com/pvthinker/pyRSW)

The configuration of JupiRSW is based on:

*Li C, Ingersoll AP, Klipfel AP, Brettle H. Modeling the stability of polygonal patterns of vortices at the poles of Jupiter as revealed by the Juno spacecraft. Proc Natl Acad Sci U S A. 2020 Sep 29;117(39):24082-24087. doi: 10.1073/pnas.2008440117. Epub 2020 Sep 8. PMID: 32900956; PMCID: PMC7533696.*

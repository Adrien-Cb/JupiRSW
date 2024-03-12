# JupiRSW

This is a rotating shallow water model for studying vortices at a planet's pole.




https://github.com/Adrien-Cb/JupiRSW/assets/130017358/60b82aa1-e99a-4bf2-a4de-1a9c1629c67a



🚧 *Work in progress* 🚧

## 🔎 Functionalities

- Save/load results as NetCDF.
- Plot animations, save them as MP4 or GIF.
- Set initial vortex latitudes and number.

This is for example the code used to generate the demo video:
```python
from model import Model
import plot

m = Model(300)
m.initialize(vort_lat=[83, 70], vort_number=[5, 1])
m.run('3 months')

plot.show_var_anim(m, 'pv', title='PV', save_as='mp4', filename='demo')
```


## 🛠 Install

### Libraries 

JupiRSW uses **NumPy**, **Xarray**, **SciPy** and **Matplotlib**. 

### ffmpeg

**ffmpeg** is required to save animations as MP4. 
There are many builds available, you can for instance get one [here](https://www.gyan.dev/ffmpeg/builds/).
Once installed, add the path to the ffmpeg executable in `config.py`:
```python
FFMPEG_PATH = R'C:\ffmpeg\bin\ffmpeg.exe'
```
Animations can still be saved as GIFs without ffmpeg.

### Optional install

I recommend using the [**cmcrameri**](https://pypi.org/project/cmcrameri/) module for pretty, perceptually uniform, colorblind-friendly colormaps. 
This install is optional.


## 📚 References and credits

The RK3 scheme implementation is taken from [PyRSW](https://github.com/pvthinker/pyRSW).

The configuration of JupiRSW is based on:

*Li C, Ingersoll AP, Klipfel AP, Brettle H. Modeling the stability of polygonal patterns of vortices at the poles of Jupiter as revealed by the Juno spacecraft. Proc Natl Acad Sci U S A. 2020 Sep 29;117(39):24082-24087. doi: 10.1073/pnas.2008440117. Epub 2020 Sep 8. PMID: 32900956; PMCID: PMC7533696.*

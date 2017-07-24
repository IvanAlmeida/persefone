PLUTO Tools
============

Assorted methods and classes to handle and visualize the output of the PLUTO MHD code.

# Simple examples 
	
#### To read the simulation output for snapshot 10 (e.g. rho.0010.dbl etc):
	
```python 
import pluto
c=pluto.Pluto(10)
```

The *p* object's attributes are now:
	
- grid: *x1,x2,x3*
- velocities: *v1,v2,v3*
- pressure *p*
- density *rho*
- number of grid cells: *n1,n2,n3*


#### Plots density field in cartesian coordinates:

```python	
c.snap()
```

or if you don't want to define an object:

```python	
pluto.Pluto(10).snap()
```

#### Plots a density field for snapshot 50 which was generated in *polar coordinates*, "regridded" with a 400x400 cartesian grid:

```python	
p=pluto.Pluto(50).pol2cart(400)
p.snap()
```
or
```python	
pluto.Pluto(50).pol2cart(400).snap()
```

#### Generate a movie of a simulation

You need to generate the image files of each frame.  For example, one quick and dirty way of doing that is to edit the *pluto.movie* method and customize what specific variable or rendering of the simulation you want to animate in the loop. 

In the directory with the image files (e.g. plot.0001.jpeg, plot.0002.jpeg etc) then run:

```shell
sh PATH/mencoder.sh 5 movie
```

This will create two movie files -- movie.avi and movie.mov (OS X compatible) -- at 5 frames/second. 

*mencoder.sh* has two dependencies: mencoder and ffmpeg.

#### More examples coming soon.

# Full torus simulation and analysis

Please see Jupyter notebook `torus_tutorial.ipynb`.

# TODO

- [ ] compute total mass in volume
- [ ] compute Mdot for a given radius
- [ ] compute energy, a.m. 


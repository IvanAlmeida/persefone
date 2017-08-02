Persefone
============

Assorted methods and classes to handle and visualize the output of the PLUTO MHD code. Adaptation in develop of the PLUTO-tools.

# Simple examples 
	
#### To read the simulation output for snapshot 10 (e.g. rho.0010.dbl etc):
	
```python 
import persefone
c=persefone.Persefone(10)
```

The *p* object's attributes are now:
	
- grid: *x1,x2,x3*
- velocities: *v1,v2,v3*
- pressure *p*
- density *rho*
- number of grid cells: *n1,n2,n3*


#### Plots density field in cartesian coordinates:

```python	
c.snapshot()
```

or if you don't want to define an object:

```python	
persefone.Persefone(10).snapshot()
```

#### Plots a density field for snapshot 50 which was generated in *polar coordinates*, "regridded" with a 400x400 cartesian grid:

```python	
p=persefone.Persefone(50).pol2cart(400)
p.snap()
```
or
```python	
persefone.Persefone(50).pol2cart(400).snapshot()
```

#### Generate a movie of a simulation

You need to generate the image files of each frame.  For example, one quick and dirty way of doing that is to edit the *pluto.movie* method and customize what specific variable or rendering of the simulation you want to animate in the loop. 

In the directory with the image files (e.g. plot.0001.jpeg, plot.0002.jpeg etc) then run:

```shell
rename 's/plot.//' *.<extension, jpeg or png>
avconv -f image2 -i %05d.png output.mp4 <it can be avi, mkv, etc>
```

This will create a movie file -- output.mp4

#### More examples coming soon.

# Full torus simulation and analysis

Please see Jupyter notebook `torus_tutorial.ipynb`.

# TODO

- [ ] compute total mass in volume
- [ ] compute Mdot for a given radius
- [ ] compute energy, a.m. 


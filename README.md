polyquad
========
Frugal generation of quadrature for polytopal domains

About
-----
polyquad is a tool that generates quadratures (or cubatures) for
polytopal domains. Important points:
- the distribution of points is shape independent
- the number of integration points depends **only** on the polynomial
  order of the integrand
- expensive operations (QR decomposition) is done only once per
  polynomial order, their results are reused from one shape to the
  other
- integration points are defined on a bounding box encapsulating the
  polytope (and can thus fall outside of the polytope), so the
  integrand should be defined there as well.


## Installation

The package is deployed on pypi, so it can be installed simply using pip

'''
pip install polyquad
'''

## Usage
For a basic usage you have at your disposal 2 functions:
- `polyquad.get_quadrature_2d`
- `polyquad.get_quadrature_3d`

These should behave nicely as long as you feed them the right data
structures. Let us explain them here

### `polyquad.get_quadrature_2d`
As the name suggests, the function should be used to get quadratures
over a 2d polygon. Let us tackle the case of a simple pentagon:

We first need to declare the coordinates of vertices

`
verts = np.array((1,-1),
                 (-1,0),
				 (0,3),
				 (2,3),
				 (3,0))
`

Then we need to specify the ordering of vertices, which is quite
straightforward in that case

`
face = np.array((0,1,2,3,4))
`

Then, getting the quadrature for a polynomial order `k` is as simple as

`
points, weights = polyquad.get_quadrature_2d(k, verts, face)
`

### `polyquad.get_quadrature_3d`
This one is for polyhedra. The only difference lies the declaration of
faces. For sake of simplicity we here give the example of a simple
pyramid (this library is intend to be used with way more complex
polyhedra)

Again we first declare the coordinates of vertices

`
verts = np.array(((0,   0, 0),
                  (1,   0, 0),
                  (1,   1, 0),
                  (0,   1, 0),
                  (.5, .5, 1)))
`

Then we need to declare each planar face. A planar face is defined by
the list of the index of vertices comprising it. For example the
square base of the pyramid contains the four first vertices defined
above, so it is defined as : `[0,1,2,3]`
This is done for each face and gathered in a list of list as follow:

`
faces = [[0,1,2,3],
         [0,1,4],
         [1,2,4],
         [2,3,4],
         [3,0,4]]
`

**Note**: if the faces of the polyhedron all have the same number of
vertices then `faces` can be casted to a `np.array` of shape (number of
faces, number of vertices per face). This slightly speeds up the
computation.

## cite this work

If you use this code for you work, you are kindly invited to cite the
associated paper
```
@article{langlois2024,
	title: {Frugal numerical integration scheme for polytopal domains},
	author: {Langlois,C. and van Putten,T. and B\'eriot,H. and Deckers,E.}
	}
```

<hr/> 

Feedback and comments can be addressed to corresponding authors
of the paper.

"""

              __                       __
   ___  ___  / /_ _____ ___ _____ ____/ /
  / _ \/ _ \/ / // / _ `/ // / _ `/ _  / 
 / .__/\___/_/\_, /\_, /\_,_/\_,_/\_,_/  
/_/          /___/  /_/                  


polyquad is a package designed to generate numerical quadratures for polytopal domains.

It is a direct implementation of the paper /Frugal numerical integration scheme for polytopal domains/,
Authors: Christophe Langlois, Thijs van Putten, Hadrien B{\'e}riot, Elke Deckers
DOI: <INSERT DOI>

keywords: polytopal domains, polygon, polyhedron, quadrature, numerical integration.

Developped and maintained by christophe and thijs
"""

from ._mapping import map_to_local_bb_2d, map_to_local_bb_3d
from ._antonietti import integrateMonomialsPolygon, integrateMonomialsPolyhedron
from ._moment_matching_3d import moment_matching as mm3d
from ._moment_matching_2d import moment_matching as mm2d
from ._main import get_quadrature_2d, get_quadrature_3d

__all__ = ['get_quadrature_2d', 'get_quadrature_3d','map_to_local_bb_2d', 'map_to_local_bb_3d', 'mm2d', 'mm3d', 'integrateMonomialsPolygon', 'integrateMonomialsPolyhedron']

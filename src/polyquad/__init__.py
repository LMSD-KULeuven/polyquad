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

from ._mapping import map_to_local_bb
from ._antonietti import integrateMonomialsPolygon, integrateMonomialsPolyhedron
from ._moment_matching import moment_matching
from ._main import get_quadrature

__all__ = ['map_to_local_bb', 'moment_matching', 'integrateMonomialsPolygon', 'integrateMonomialsPolyhedron', 'get_quadrature']

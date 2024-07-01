"""
              __                       __
   ___  ___  / /_ _____ ___ _____ ____/ /
  / _ \/ _ \/ / // / _ `/ // / _ `/ _  / 
 / .__/\___/_/\_, /\_, /\_,_/\_,_/\_,_/  
/_/          /___/  /_/                  


               ___                                         __     
              /\_ \                                       /\ \    
 _____     ___\//\ \    __  __     __   __  __     __     \_\ \   
/\ '__`\  / __`\\ \ \  /\ \/\ \  /'__`\/\ \/\ \  /'__`\   /'_` \  
\ \ \L\ \/\ \L\ \\_\ \_\ \ \_\ \/\ \L\ \ \ \_\ \/\ \L\.\_/\ \L\ \ 
 \ \ ,__/\ \____//\____\\/`____ \ \___, \ \____/\ \__/.\_\ \___,_\
  \ \ \/  \/___/ \/____/ `/___/> \/___/\ \/___/  \/__/\/_/\/__,_ /
   \ \_\                    /\___/    \ \_\                       
    \/_/                    \/__/      \/_/                       



                |                                  | 
  __ \    _ \   |  |   |   _` |  |   |   _` |   _` | 
  |   |  (   |  |  |   |  (   |  |   |  (   |  (   | 
  .__/  \___/  _| \__, | \__, | \__,_| \__,_| \__,_| 
 _|               ____/      _|                      



              |                                | 
   _ \   _ \  |  |  |   _` |  |  |   _` |   _` | 
  .__/ \___/ _| \_, | \__, | \_,_| \__,_| \__,_| 
 _|             ___/      _|                     

polyquad is a python package designed to generated quadratures for polytopal domains.

It is a direct implementation of the paper /Frugal numerical integration scheme for polytopal domains/,
Authors: C.Langlois, T. van Putten, H. B{\'e}riot, E.Deckers
DOI: <INSERT DOI>

keywords: polytopal domains, polygon, polyhedron, quadrature, numerical integration.

Developped by christophe and thijs
"""

from ._mapping import map_to_local_bb
from ._antonietti import integrate_monomials_polyhedron
from ._moment_matching import moment_matching
from ._main import get_quadrature

__all__ = ['map_to_local_bb', 'moment_matching', 'integrate_monomials_polyhedron', 'get_quadrature']

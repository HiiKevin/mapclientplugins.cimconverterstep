
"""
MAP Client Plugin - Generated from MAP Client v0.16.0
"""

__version__ = '0.1.0'
__author__ = 'Xxxx Yyyyy'
__stepname__ = 'CIM Converter'
__location__ = ''

# import class that derives itself from the step mountpoint.
from mapclientplugins.cimconverterstep import step

# Import the resource file when the module is loaded,
# this enables the framework to use the step icon.
from . import resources_rc
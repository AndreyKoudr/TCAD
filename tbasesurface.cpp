#include "tbasesurface.h"

/*
  UNDEFINED = -1,
  INACTIVE,
  WING,
  CYLINDER,
  CYLINDER_TOP,
  CYLINDER_BOTTOM,
  CYLINDER_FLANK0,
  CYLINDER_FLANK1,
  SPHERE,
  ELLIPSOID,
  BOX,
  SHAFT,
  CASING,

  FACETYPE_TOTAL
*/

std::string tcad::facetypeNames[FACETYPE_TOTAL] = 
{
  "undefined",
  "inactive",
  "wing",
  "cylinder",
  "cylinder_top",
  "cylinder_bottom",
  "cylinder_flank0",
  "cylinder_flank1",
  "sphere",
  "ellipsoid",
  "box",
  "shaft",
  "casing"
};


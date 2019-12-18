#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/imagesSetsUtils/SetFromImage.h>
#include <DGtal/images/ImageSelector.h>
#include <DGtal/images/ConstImageAdapter.h>
#include <DGtal/images/RigidTransformation2D.h>
#include <DGtal/io/boards/Board2D.h>
#include <DGtal/io/Color.h>
#include <DGtal/io/colormaps/ColorBrightnessColorMap.h>
#include <DGtal/io/readers/PGMReader.h>
#include <DGtal/io/writers/GenericWriter.h>
#include <DGtal/topology/SurfelAdjacency.h>
#include <DGtal/topology/helpers/Surfaces.h>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <cmath>
#include <iostream>
#include <string>

extern const std::string PATH = "../TP3_images/binary/";
extern const std::string CLE = "Cle_";
extern const std::string PINCE = "Pince_";
extern const std::string CUTTER = "Cutter_";
extern const std::string FORMAT = ".pgm";

#include "Find2DAffineTransform.cpp"
#include "Component.cpp"

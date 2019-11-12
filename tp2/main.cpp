#include <DGtal/base/Common.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/ImageSelector.h>
#include "DGtal/io/readers/PGMReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include <DGtal/images/imagesSetsUtils/SetFromImage.h>
#include <DGtal/io/boards/Board2D.h>
#include <DGtal/io/colormaps/ColorBrightnessColorMap.h>
#include <DGtal/topology/SurfelAdjacency.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include "DGtal/io/Color.h"

using namespace std;
using namespace DGtal;
using namespace Z2i;


template <class T>
size_t getPerimeter(T &object) {
    // make a Kovalevsky-Khalimsky space
    KSpace t_KSpace;
    t_KSpace.init(object.domain().lowerBound() - Point(2,2), object.domain().upperBound() + Point(2,2), true);
    // set an adjacency (4-connectivity)
    SurfelAdjacency<2> sAdj(true);

    // search for one boundary element
    SCell bel = Surfaces<KSpace>::findABel(t_KSpace, object.pointSet(), 10000);

    // boundary tracking
    std::vector<Z2i::Point> t_BoundaryPoints;
    Surfaces<Z2i::KSpace>::track2DBoundaryPoints(t_BoundaryPoints, t_KSpace, sAdj, object.pointSet(), bel);

    return t_BoundaryPoints.size();

}

template <class T>
size_t getArea(T &object) {

    return object.size();
}

template <class T>
Curve getBoundary(T & object)
{
    // make a Kovalevsky-Khalimsky space
    KSpace t_KSpace;
    t_KSpace.init(object.domain().lowerBound() - Point(2,2), object.domain().upperBound() + Point(2,2), true);
    // set an adjacency (4-connectivity)
    SurfelAdjacency<2> sAdj(true);

    // search for one boundary element
    SCell bel = Surfaces<KSpace>::findABel(t_KSpace, object.pointSet(), 10000);

    // boundary tracking
    std::vector<Z2i::Point> t_BoundaryPoints;
    Surfaces<Z2i::KSpace>::track2DBoundaryPoints(t_BoundaryPoints, t_KSpace, sAdj, object.pointSet(), bel);

    // obtain a curve
    Curve boundaryCurve;
    boundaryCurve.initFromVector(t_BoundaryPoints);

    return boundaryCurve;
}


int main(int argc, char** argv)
{
    typedef ImageSelector<Domain, unsigned char>::Type Image;
    typedef Domain::ConstIterator DomainConstIterator;
    typedef DigitalSetSelector<Domain, BIG_DS+HIGH_BEL_DS>::Type DigitalSet;
    typedef Object<DT4_8, DigitalSet> ObjectType4_8;
    typedef Object<DT8_4, DigitalSet> ObjectType8_4;

    // read an image
    Image image = PGMReader<Image>::importPGM (argv[1]);

    // make a digital set from the image
    Z2i::DigitalSet set2d(image.domain());                                // Create a digital set of proper size
    SetFromImage<Z2i::DigitalSet>::append<Image>(set2d, image, 1, 255);    // Populate a digital set from the input image

    // create a digital object from the digital set
    std::vector<ObjectType4_8> objects4_8;     // All connected components are going to be stored in it
    std::back_insert_iterator<std::vector<ObjectType4_8>> inserter4_8(objects4_8);

    // create a digital object from the digital set
    std::vector<ObjectType8_4> objects8_4;     // All connected components are going to be stored in it
    std::back_insert_iterator<std::vector<ObjectType8_4>> inserter8_4(objects8_4);

    // obtain the connected components with a chosen adjacency pair
    ObjectType4_8 bdiamond4_8(dt4_8, set2d);      // (4,8) adjacency
    bdiamond4_8.writeComponents(inserter4_8);

    // obtain the connected components with a chosen adjacency pair
    ObjectType8_4 bdiamond8_4(dt8_4, set2d);      // (8,4) adjacency
    bdiamond8_4.writeComponents(inserter8_4);

    // for each grain, make the boundary tracking, calculate the area, perimeter, etc....

    int index = 0;
    float areaSum = 0;
    float count = objects4_8.size();
    float areaMean;
    for(auto object4_8 : objects4_8) {
      try {
        areaSum += getArea(object4_8);
      } catch(std::exception e) {
        std::cout << "at [" << index << "] : " << e.what() << std::endl;
        --count;
      }
      ++index;
    }

    areaMean = areaSum / count;

    std::cout << "area mean : " << areaMean << std::endl;

    //stock indexes of valid grains
    std::vector<int> indexes;
    index = 0;
    for(auto object4_8 : objects4_8) {
      try {
        float area = getArea(object4_8);
        if(area >= areaMean) {
          indexes.push_back(index);
        }
      } catch(std::exception e) {
        std::cout << "at [" << index << "] : " << e.what() << std::endl;
      }
      ++index;
    }

    //for each valid grain
    for(auto i : indexes) {
      try {
        Z2i::Curve c = getBoundary(objects4_8[i]);

        // make a pdf file of the object
        Board2D aBoard;
        std::stringstream filename;
        filename << "export/boundary_" << i << ".pdf";
        aBoard << objects4_8[i];
        aBoard << c;
        aBoard.saveCairo(filename.str().c_str(), Board2D::CairoPDF);
      } catch(std::exception e) {
        std::cout << "at [" << i << "] : " << e.what() << std::endl;
      }
    }

		std::cout << "Grains count in (4,8) adjacency : " << objects4_8.size() << std::endl;
		std::cout << "Grains count in (8,4) adjacency : " << objects8_4.size() << std::endl;
		std::cout << "well-composed : " << ((objects4_8.size() != objects8_4.size()) ? "false" : "true") << std::endl;

    return 0;
}

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

#include <Eigen/Dense>
#include <Eigen/Geometry>

using namespace std;
using namespace DGtal;
using namespace Z2i;

const string PATH = "../TP3_images/binary/";
const string CLE = "Cle_";
const string PINCE = "Pince_";
const string CUTTER = "Cutter_";
const string FORMAT = ".pgm";

typedef ImageSelector<Domain, unsigned char >::Type TypeImage;
typedef Domain::ConstIterator DomainConstIterator;
typedef DigitalSetSelector< Domain, BIG_DS+HIGH_BEL_DS >::Type DigitalSet;
typedef Object<DT4_8, DigitalSet> ObjectType4_8;

Eigen::Affine2d Find2DAffineTransform(const Eigen::Matrix2Xd &_in, const Eigen::Matrix2Xd &_out);

std::string toString(const Eigen::MatrixXd& mat){
    std::stringstream ss;
    ss << mat;
    return ss.str();
}

template <class T>
Curve getBoundary(T & object)
{
    // make a Kovalevsky-Khalimsky space
    KSpace t_KSpace;
    t_KSpace.init( object.domain().lowerBound() - Point(2,2), object.domain().upperBound() + Point(2,2), true );
    // set an adjacency (4-connectivity)
    SurfelAdjacency<2> sAdj( true );

    // search for one boundary element
    SCell bel = Surfaces<KSpace>::findABel(t_KSpace, object.pointSet(), 10000);

    // boundary tracking
    std::vector<Z2i::Point> t_BoundaryPoints;
    Surfaces<Z2i::KSpace>::track2DBoundaryPoints( t_BoundaryPoints, t_KSpace, sAdj, object.pointSet(), bel );

    // obtain a curve
    Curve boundaryCurve;
    boundaryCurve.initFromVector( t_BoundaryPoints );

    return boundaryCurve;
}

template<class T>
Z2i::Point getCenterOfMass(T object) {
    Z2i::Point p(0, 0);
    for(auto it = object.begin(), itEnd = object.end(); it != itEnd; ++it) {
        std::cout << *it << std::endl;
        p += *it;
    }

    return p / int(object.size());
}

class Component {
    public:
    Component() = delete;
    Component(std::string filePath) {
      // read an image
      TypeImage image = PGMReader<TypeImage>::importPGM (filePath);

      // make a digital set from the image
      Z2i::DigitalSet set2d (image.domain());                                 // Create a digital set of proper size
      SetFromImage<Z2i::DigitalSet>::append<TypeImage>(set2d, image, 1, 255);     //populate a digital set from the input image

      // create a digital object from the digital set
      std::back_insert_iterator< std::vector< ObjectType4_8 > > inserter( objects );

      // obtain the connected components with a chosen adjacency pair
      ObjectType4_8 bdiamond(dt4_8, set2d);      // (4,8) adjacency
      bdiamond.writeComponents(inserter);

      if(objects.size() == 0) {
        throw new std::out_of_range("No component found !");
      }

      std::cout << objects[0].size() << std::endl;
    }

    const Z2i::Point getCenterOfMass() const {
      Z2i::Point p(0, 0);
      for(auto it = objects[0].begin(), itEnd = objects[0].end(); it != itEnd; ++it) {
          p += *it;
      }

      return p / int(objects[0].size());
    }


    const Eigen::Matrix2Xd getMatrixOfPoints() const {
      Eigen::Matrix2Xd m(2, objects[0].size());
      int index = 0;
      for(auto it = objects[0].begin(), itEnd = objects[0].end(); it != itEnd; ++it) {
          m(0, index) = double((*it)[0]);
          m(1, index) = double((*it)[1]);
          ++index;
      }
      return m;
    }

    std::vector< ObjectType4_8 > objects; // All connected components are going to be stored in it
};


int main(int argc, char** argv)
{

    Component cle0(PATH + CLE + "00" + FORMAT);
    Component cle1(PATH + CLE + "01" + FORMAT);

    Z2i::Point centerOfMass0 = cle0.getCenterOfMass();
    Z2i::Point centerOfMass1 = cle1.getCenterOfMass();

    Z2i::Point translation01 = centerOfMass1 - centerOfMass0;

    std::cout << translation01 << std::endl;

    // for each grain, make the boundary tracking, calculate the area, perimeter, etc....
    //std::cout << "Get the object boundary." << endl;
    //Z2i::Curve c = getBoundary(objects[0]);

    // make a pdf file
    Board2D aBoard;
    aBoard << cle0.objects[0];
    aBoard.setFillColor(Color::Red);
    aBoard.drawCircle(centerOfMass0[0], centerOfMass0[1], 10);
    aBoard.setFillColor(Color::Cyan);
    aBoard.drawCircle(centerOfMass1[0], centerOfMass1[1], 5);
    aBoard.saveCairo("boundaryCurve.pdf",Board2D::CairoPDF);

    Eigen::Matrix3d affine = Find2DAffineTransform(cle0.getMatrixOfPoints(), cle1.getMatrixOfPoints()).matrix();

    std::cout << affine << std::endl;

    return 0;
}


// This code is released in public domain

// Given two sets of 3D points, find the rotation + translation + scale
// which best maps the first set to the second.
// Source: http://en.wikipedia.org/wiki/Kabsch_algorithm

//https://github.com/oleg-alexandrov/projects/blob/master/eigen/Kabsch.cpp

// The input 2D points are stored as columns.
Eigen::Affine2d Find2DAffineTransform(const Eigen::Matrix2Xd &_in, const Eigen::Matrix2Xd &_out) {

  // Default output
  Eigen::Affine2d A;
  A.linear() = Eigen::Matrix2d::Identity(2, 2);
  A.translation() = Eigen::Vector2d::Zero();

  Eigen::Matrix2Xd in = _in;
  Eigen::Matrix2Xd out = _out;

  if(in.cols() > out.cols()) { 
    in.resize(2, out.cols());
  } else {
    out.resize(2, in.cols());
  }

  std::cout << in.cols() << " " << out.cols() << std::endl;

  // First find the scale, by finding the ratio of sums of some distances,
  // then bring the datasets to the same scale.
  double dist_in = 0, dist_out = 0;
  for (int col = 0; col < in.cols()-1; col++) {
    dist_in  += (in.col(col+1) - in.col(col)).norm();
    dist_out += (out.col(col+1) - out.col(col)).norm();
  }
  if (dist_in <= 0 || dist_out <= 0)
    return A;
  double scale = dist_out/dist_in;
  out /= scale;

  // Find the centroids then shift to the origin
  Eigen::Vector2d in_ctr = Eigen::Vector2d::Zero();
  Eigen::Vector2d out_ctr = Eigen::Vector2d::Zero();
  for (int col = 0; col < in.cols(); col++) {
    in_ctr  += in.col(col);
    out_ctr += out.col(col);
  }
  in_ctr /= in.cols();
  out_ctr /= out.cols();
  for (int col = 0; col < in.cols(); col++) {
    in.col(col)  -= in_ctr;
    out.col(col) -= out_ctr;
  }

  // SVD
  Eigen::MatrixXd Cov = in * out.transpose();
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

  // Find the rotation
  double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
  if (d > 0)
    d = 1.0;
  else
    d = -1.0;
  Eigen::Matrix2d I = Eigen::Matrix2d::Identity(2, 2);
  I(1, 1) = d;
  Eigen::Matrix2d R = svd.matrixV() * I * svd.matrixU().transpose();

  // The final transform
  A.linear() = scale * R;
  A.translation() = scale*(out_ctr - R*in_ctr);

  return A;
}

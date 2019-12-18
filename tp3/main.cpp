#include "./main.hpp"

using namespace std;
using namespace DGtal;
using namespace Z2i;

const string PATH = "../TP3_images/binary/";
const string CLE = "Cle_";
const string PINCE = "Pince_";
const string CUTTER = "Cutter_";
const string FORMAT = ".pgm";

typedef ImageSelector<Domain, unsigned char>::Type TypeImage;
typedef Domain::ConstIterator DomainConstIterator;
typedef DigitalSetSelector<Domain, BIG_DS + HIGH_BEL_DS>::Type DigitalSet;
typedef Object<DT4_8, DigitalSet> ObjectType4_8;

std::string toString(const Eigen::MatrixXd &mat)
{
  std::stringstream ss;
  ss << mat;
  return ss.str();
}

template <class T>
Curve getBoundary(T &object)
{
  // make a Kovalevsky-Khalimsky space
  KSpace t_KSpace;
  t_KSpace.init(object.domain().lowerBound() - Point(2, 2), object.domain().upperBound() + Point(2, 2), true);
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

int main(int argc, char **argv)
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
  aBoard.saveCairo("boundaryCurve.pdf", Board2D::CairoPDF);

  Eigen::Matrix3d affine = Find2DAffineTransform(cle0.getMatrixOfPoints(), cle1.getMatrixOfPoints()).matrix();

  std::cout << affine << std::endl;

  return 0;
}

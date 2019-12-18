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

string toString(const Eigen::MatrixXd &mat)
{
  stringstream ss;
  ss << mat;
  return ss.str();
}

const Z2i::Point getTranslationBetweenComponents(Component &c1, Component &c2)
{
  Z2i::Point centerOfMass0 = c1.getCenterOfMass();
  Z2i::Point centerOfMass1 = c2.getCenterOfMass();
  Z2i::Point translation = centerOfMass1 - centerOfMass0;
  return translation;
}

const Eigen::Matrix2d getRotationBetweenComponents(Component &c1, Component &c2)
{
  Eigen::Matrix3d affine = Find2DAffineTransform(c1.getMatrixOfPoints(), c2.getMatrixOfPoints()).matrix();
  Eigen::Matrix2d rotation = affine.topLeftCorner(2, 2);
  return rotation;
}

int main(int argc, char **argv)
{

  Component cle0(PATH + CLE + "00" + FORMAT);
  Component cle1(PATH + CLE + "01" + FORMAT);

  cout << "Translation between cle0 and cle1 is " << getTranslationBetweenComponents(cle0, cle1) << endl;

  // make a pdf file

  // Step 3::1
  // Display the translation vector

  /*
  Board2D aBoard;
  aBoard << cle0.largestObject;
  aBoard.setFillColor(Color::Red);
  Z2i::Point centerOfMass0 = cle0.getCenterOfMass();
  aBoard.drawCircle(centerOfMass0[0], centerOfMass0[1], 10);
  Z2i::Point centerOfMass1 = cle1.getCenterOfMass();
  aBoard.setFillColor(Color::Cyan);
  aBoard.drawCircle(centerOfMass1[0], centerOfMass1[1], 5);
  aBoard.saveCairo("boundaryCurve.pdf", Board2D::CairoPDF);

  */

  // Step 3::2
  // Display the rotation matrix
  cout << "Rotation is " << getRotationBetweenComponents(cle0, cle1) << endl;

  return 0;
}

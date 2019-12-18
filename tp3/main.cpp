#include "./main.hpp"

using namespace std;
using namespace DGtal;
using namespace functors;
using namespace Z2i;

typedef ImageSelector<Domain, unsigned char>::Type TypeImage;
typedef Domain::ConstIterator DomainConstIterator;
typedef DigitalSetSelector<Domain, BIG_DS + HIGH_BEL_DS>::Type DigitalSet;
typedef Object<DT4_8, DigitalSet> ObjectType4_8;

// Step 4
// Define types for rigid transformation

typedef ForwardRigidTransformation2D<Space> ForwardTrans;
typedef BackwardRigidTransformation2D<Space> BackwardTrans;
typedef DomainRigidTransformation2D<Domain, ForwardTrans> DomainTransformer;
typedef ConstImageAdapter<TypeImage, Domain, BackwardTrans, TypeImage::Value, Identity> ImageBackwardAdapter;
typedef DomainTransformer::Bounds Bounds;

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

  // make a pdf file

  // Step 3::1
  // Display the translation vector

  Z2i::Point translation = getTranslationBetweenComponents(cle0, cle1);
  cout << "Translation between cle0 and cle1 is " << translation << endl;

  Board2D aBoard;
  aBoard << cle0.largestObject;
  aBoard.setFillColor(Color::Red);
  Z2i::Point centerOfMass0 = cle0.getCenterOfMass();
  aBoard.drawCircle(centerOfMass0[0], centerOfMass0[1], 10);
  Z2i::Point centerOfMass1 = cle1.getCenterOfMass();
  aBoard.setFillColor(Color::Cyan);
  aBoard.drawCircle(centerOfMass1[0], centerOfMass1[1], 5);
  aBoard.saveCairo("boundaryCurve.pdf", Board2D::CairoPDF);

  // Step 3::2
  // Display the rotation matrix
  Eigen::Matrix2d rotation = getRotationBetweenComponents(cle0, cle1);
  cout << "Rotation between cle0 and cle1 is " << endl;
  cout << rotation << endl;
  float angleInRadian = acos(rotation(0, 0));
  cout << "Angle is " << angleInRadian << endl;

  // Step 4
  // Create transformation by providing information about origin, angle in radians, and translation
  cle1.createImageBackwardAdapter(RealPoint(0, 0), -angleInRadian, translation);

  return 0;
}

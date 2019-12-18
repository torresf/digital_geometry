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

// Step 5
// Distance transform
typedef IntervalForegroundPredicate<TypeImage> Binarizer;
typedef DistanceTransformation<Space, Binarizer, L2Metric> DTL2;
typedef HueShadeColorMap<DTL2::Value, 2> HueTwice;

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

// Inspired by this article on stackoverflow
// https://stackoverflow.com/questions/46637401/hausdorff-distance-object-detection

const int squaredEuclideanDistance(const Z2i::Point &p1, const Z2i::Point &p2)
{
  Z2i::Point d = p1 - p2;
  return (d[0] * d[0]) + (d[1] * d[1]);
}

const int maxPointDistance(const Component &c1, const Component &c2)
{
  // Find the maximum distance
  int maxDistance = 0;
  for (Z2i::Point point : c1.largestObject)
  {
    // Find the minimum distance
    int minDistance = numeric_limits<int>::max();

    for (Z2i::Point tempPoint : c2.largestObject)
    {
      int distance = squaredEuclideanDistance(point, tempPoint);
      if (distance < minDistance)
      {
        minDistance = distance;
      }
      if (distance == 0)
      {
        break;
      }
    }
    maxDistance += minDistance;
  }

  return maxDistance;
}

// Step 5::1
const float hausdorffDistance(const Component &c1, const Component &c2)
{
  // Hausdorff distance between two components
  // Defined by dH(S1,S′2)=max{supx∈S1{infy∈S′2d(x,y)}, supy∈S′2{infx∈S1d(x,y)}}

  float maxDistc1 = maxPointDistance(c1, c2);
  float maxDistc2 = maxPointDistance(c2, c1);

  return sqrt(max(maxDistc1, maxDistc2));
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
  const string reconstructedFilename = "cle1_transformed";
  cle1.createImageBackwardAdapter(RealPoint(0, 0), -angleInRadian, translation, "");

  // Step 5
  // Compare the transformed image and the reference
  Component reconstructedCle0("../TP3_images/created/" + reconstructedFilename + FORMAT);
  //float distance = hausdorffDistance(cle0, cle1);
  //cout << "Hausdorff distance is " << distance << endl;

  DTL2 dt = reconstructedCle0.getDTL2();

  return 0;
}

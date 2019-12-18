#include "Component.hpp"

using namespace std;
using namespace DGtal;
using namespace Z2i;

typedef ImageSelector<Domain, unsigned char>::Type TypeImage;
typedef Domain::ConstIterator DomainConstIterator;
typedef DigitalSetSelector<Domain, BIG_DS + HIGH_BEL_DS>::Type DigitalSet;
typedef Object<DT4_8, DigitalSet> ObjectType4_8;

class Component
{
public:
  Component() = delete;
  Component(string filePath)
  {

    // STEP 2::1
    // Use PGMReader to read each input image.
    // read an image
    TypeImage image = PGMReader<TypeImage>::importPGM(filePath);

    // STEP 2::2
    // Convert a binary image to a "digital set" following the instruction of "DigitalSet from threholded image".
    // make a digital set from the image
    Z2i::DigitalSet set2d(image.domain());                                  // Create a digital set of proper size
    SetFromImage<Z2i::DigitalSet>::append<TypeImage>(set2d, image, 1, 255); //populate a digital set from the input image

    // STEP 2::3
    // Construct a "digital object" from a "digital set" with a choice of a good adjacency pair.
    // create a digital object from the digital set
    back_insert_iterator<vector<ObjectType4_8>> inserter(this->objects);

    // STEP 2::4
    // Compute connected components of a "digital object" using WriteComponents just in case. If there are several components, take the largest one.
    // obtain the connected components with a chosen adjacency pair
    ObjectType4_8 bdiamond(dt4_8, set2d); // (4,8) adjacency
    bdiamond.writeComponents(inserter);

    if (this->objects.size() == 0)
    {
      throw new out_of_range("No component found !");
    }

    cout << "Amount of this->objects is " << this->objects.size() << endl;

    // STEP 2::4
    // If there are several components, take the largest one.
    this->largestObject = this->getLargestObject();
    cout << "Largest object pixels amount is " << this->largestObject.size() << endl;

  }

  const ObjectType4_8 getLargestObject() const
  {
    ObjectType4_8 largestObject;
    if (this->objects.size() > 1) {
      int maxArea = 0;
      for (auto o : this->objects) {
        if (o.size() > maxArea) {
          maxArea = o.size();
          largestObject = o;
        }
      }
    } else {
      largestObject = this->objects[0];
    }
    return largestObject;
  }

  const Z2i::Point getCenterOfMass() const
  {
    Z2i::Point p(0, 0);
    for (auto it = this->largestObject.begin(), itEnd = this->largestObject.end(); it != itEnd; ++it)
    {
      p += *it;
    }

    return p / int(this->largestObject.size());
  }

  const Eigen::Matrix2Xd getMatrixOfPoints() const
  {
    Eigen::Matrix2Xd m(2, this->largestObject.size());
    int index = 0;
    for (auto it = this->largestObject.begin(), itEnd = this->largestObject.end(); it != itEnd; ++it)
    {
      m(0, index) = double((*it)[0]);
      m(1, index) = double((*it)[1]);
      ++index;
    }
    return m;
  }

  vector<ObjectType4_8> objects; // All connected components are going to be stored in it
  ObjectType4_8 largestObject;
};

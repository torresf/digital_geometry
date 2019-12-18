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

class Component
{
public:
  Component() = delete;
  Component(string filePath)

      // STEP 2::1
      // Use PGMReader to read each input image.
      // read an image
      : image(PGMReader<TypeImage>::importPGM(filePath))
  {

    // STEP 2::2
    // Convert a binary image to a "digital set" following the instruction of "DigitalSet from threholded image".
    // make a digital set from the image
    // Create a digital set of proper size
    Z2i::DigitalSet set2d(image.domain());
    //populate a digital set from the input image
    SetFromImage<Z2i::DigitalSet>::append<TypeImage>(set2d, image, 1, 255);

    // STEP 2::3
    // Construct a "digital object" from a "digital set" with a choice of a good adjacency pair.
    // create a digital object from the digital set
    back_insert_iterator<vector<ObjectType4_8>> inserter(objects);

    // STEP 2::4
    // Compute connected components of a "digital object" using WriteComponents just in case. If there are several components, take the largest one.
    // obtain the connected components with a chosen adjacency pair
    ObjectType4_8 bdiamond(dt4_8, set2d); // (4,8) adjacency
    bdiamond.writeComponents(inserter);

    if (objects.size() == 0)
    {
      throw new out_of_range("No component found !");
    }

    cout << "Amount of objects is " << objects.size() << endl;

    // STEP 2::4
    // If there are several components, take the largest one.
    this->largestObject = this->getLargestObject();
    cout << "Largest object pixels amount is " << this->largestObject.size() << endl;
  }

  const ObjectType4_8 getLargestObject() const
  {
    ObjectType4_8 largestObject;
    if (this->objects.size() > 1)
    {
      int maxArea = 0;
      for (auto o : this->objects)
      {
        if (o.size() > maxArea)
        {
          maxArea = o.size();
          largestObject = o;
        }
      }
    }
    else
    {
      largestObject = this->objects[0];
    }
    return largestObject;
  }

  const Z2i::Point getCenterOfMass() const
  {
    Z2i::Point p(0, 0);
    for (auto it = largestObject.begin(), itEnd = largestObject.end(); it != itEnd; ++it)
    {
      p += *it;
    }

    return p / int(largestObject.size());
  }

  const Eigen::Matrix2Xd getMatrixOfPoints() const
  {
    Eigen::Matrix2Xd m(2, largestObject.size());
    int index = 0;
    for (auto it = largestObject.begin(), itEnd = largestObject.end(); it != itEnd; ++it)
    {
      m(0, index) = double((*it)[0]);
      m(1, index) = double((*it)[1]);
      ++index;
    }
    return m;
  }

  const ImageBackwardAdapter createImageBackwardAdapter(const RealPoint &origin,
                                                        float angleInRadian,
                                                        const RealVector &translation,
                                                        const string imageName = "backward_transform") const
  {
    cout << origin << endl;
    cout << angleInRadian << endl;
    cout << translation << endl;
    const Identity idD;
    const ForwardTrans forwardTrans(origin, angleInRadian, translation);
    const BackwardTrans backwardTrans(origin, angleInRadian, translation);
    const DomainTransformer domainTransformer(forwardTrans);
    const Bounds bounds = domainTransformer(image.domain());
    const Domain transformedDomain(bounds.first, bounds.second);
    const ImageBackwardAdapter backwardImageAdapter(image, transformedDomain, backwardTrans, idD);
    backwardImageAdapter >> "../TP3_images/created/" + imageName + FORMAT;
    cout << "Image created" << endl;
    return backwardImageAdapter;
  }

  vector<ObjectType4_8> objects; // All connected components are going to be stored in it
  ObjectType4_8 largestObject;
  TypeImage image;
};

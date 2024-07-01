#include <iostream>
#include <fstream> // Use to read/write to files
#include <sstream>
#include <string>
#include <vector>

class Point
{
public:
    int x;
    int y;

    Point(int x, int y) : x(x), y(y) {}

    bool operator==(const Point& otherPoint) const
    {
        return x == otherPoint.x && y == otherPoint.y;
    }
};

/**
 * @brief Represents a geometric shape defined by an outline and optional cutout regions.
 *
 * The Shape class provides methods for:
 * - Reading shape data from a file and constructing a Shape object.
 * - Checking if a point is inside the shape (including on the outline or cutout edges).
 * - Managing and modifying the outline and cutout vertices.
 *
 * The outline vertices define the boundary of the shape, and the cutouts are internal regions
 * in which points are not considered to be inside the shape.
 *
 * @note Uses std::move to construct/update, ensure temporary variables are passed.
 * @note Point-in-polygon algorithms (Ray Casting) are used for geometric checks.
 */
class Shape
{
private:
    std::vector<Point> outlineVertices;
    std::vector<std::vector<Point>> cutouts;

public:
    /**
     * Constructs a Shape object with given outline vertices and optional cutouts.
     *
     * @param outlineVerts The vertices defining the outline of the shape.
     * @param cutouts The internal cutout regions inside the shape.
     */
    Shape(std::vector<Point>&& outlineVerts, std::vector<std::vector<Point>>&& cutouts) : outlineVertices(std::move(outlineVerts)), cutouts(std::move(cutouts)) {}

    /**
     * Creates a Shape object based on OUTLINE and CUT data from a text file.
     *
     * @param fileName The name of the text file (must include .txt).
     * @return The Shape object.
     * @throws std::runtime_error if file failed to open.
     */
    static Shape readShapeFromFile(const std::string& fileName)
    {
        std::ifstream textFile(fileName);
        if (!textFile.is_open())
        {
            std::cerr << "Error opening " << fileName << std::endl;
            throw std::runtime_error("File cannot be opened");
        }

        // Prepare outline and cutout vectors
        std::vector<Point> shapeOutline;
        std::vector<std::vector<Point>> shapeCutouts;

        std::string line;
        while (std::getline(textFile, line)) // Refactor to read pairs at once (iss >> data1 >> data2)
        {
            std::istringstream iss(line);
            std::string shapeType;
            iss >> shapeType;

            // If current line involves keyword, read following number (number of vertices)
            // Loop through the next (number of vertices) lines to fetch vertex data
            if (shapeType == "OUTLINE")
            {
                int numVertices;
                iss >> numVertices;
                shapeOutline.reserve(numVertices); // Reserve space for vertices to reduce reallocations

                for (int i = 0; i < numVertices; ++i)
                {
                    std::getline(textFile, line);
                    std::istringstream vertexStream(line);
                    int x, y;
                    vertexStream >> x >> y;
                    shapeOutline.emplace_back(x, y); // Use emplace to construct Point object directly within vector
                }
            }
            else if (shapeType == "CUT")
            {
                int numVertices;
                iss >> numVertices;

                std::vector<Point> cutoutVertices;
                cutoutVertices.reserve(numVertices);

                for (int i = 0; i < numVertices; ++i)
                {
                    std::getline(textFile, line);
                    std::istringstream vertexStream(line);
                    int x, y;
                    vertexStream >> x >> y;
                    cutoutVertices.emplace_back(x, y);
                }

                shapeCutouts.push_back(cutoutVertices);
            }
        }

        textFile.close(); // Close the textfile

        return Shape(std::move(shapeOutline), std::move(shapeCutouts)); // Use std::move to transfer ownership of vectors (no copying in member assignment)
    }

    /**
     * Performs point-in-polygon checks using object's outline and cutout vertices.
     * A point on the edge of an outline or cutout is considered to be inside the shape.
     *
     * @param point The point (Point object) in question.
     * @return true: point is inside shape | false: point is outside shape.
     */
    bool isPointInside(const Point& point)
    {
        // First, check if the point is exactly on any of the vertices

        // Check outline vertices
        for (const auto& vertex : outlineVertices)
        {
            if (point == vertex) return true;
        }

        // Check cutout vertices
        for (const auto& cutout : cutouts)
        {
            for (const auto& vertex : cutout)
            {
                if (point == vertex) return true;
            }
        }

        // If point isn't on a vertex, carry out more complex checks
        return (isPointInsideOutline(point) && !isPointInsideCutout(point)) || (isPointOnOutlineEdge(point) || isPointOnCutoutEdge(point));
    }

    /**
     * Adds a cutout region to the shape.
     *
     * @param cutoutVerts The vertices defining the cutout region.
     */
    void addCutout(std::vector<Point>&& cutoutVerts)
    {
        cutouts.push_back(std::move(cutoutVerts));
    }

    // Getters
    const std::vector<Point> getOutlineVertices() const
    {
        return outlineVertices;
    }

    const std::vector<std::vector<Point>>& getCutouts() const
    {
        return cutouts;
    }

    // Setters
    void setOutlineVertices(std::vector<Point>&& newOutlineVertices)
    {
        outlineVertices = std::move(newOutlineVertices);
    }

    void setCutouts(std::vector<std::vector<Point>>&& newCutouts)
    {
        cutouts = std::move(newCutouts);
    }

private:
    /**
     * Evaluates whether a given point is inside the shape outline.
     *
     * @param point The point (Point object) in question.
     * @return true: the point is inside the outline | false: the point is not inside the outline.
     */
    bool isPointInsideOutline(const Point& point)
    {
        size_t numVertices = outlineVertices.size();
        bool isInside = false;
        
        // Ray Casting Algorithm
        for (size_t i = 0, j = numVertices - 1; i < numVertices; j = i++)
        {
            if (((outlineVertices[i].y > point.y) != (outlineVertices[j].y > point.y)) 
                && (point.x < (outlineVertices[j].x - outlineVertices[i].x) * (point.y - outlineVertices[i].y) / (outlineVertices[j].y - outlineVertices[i].y) + outlineVertices[i].x))
            {
                isInside = !isInside;
            }
        }
        
        return isInside;
    }

    /**
     * Evaluates whether a given point is on an edge of the shape outline.
     *
     * @param point The point (Point object) in question.
     * @param epsilon Tolerance threshold (to account for numerical precision issues). Defaults to 1e-6.
     * @return true: the point is on the edge | false: the point is not on the edge.
     */
    bool isPointOnOutlineEdge(const Point& point, float epsilon = 1e-6)
    {
        size_t numVertices = outlineVertices.size();
        float xIntersection;

        for (size_t i = 0, j = numVertices - 1; i < numVertices; j = i++)
        {
            // Check for horizontal edge (avoid division by zero)
            if (outlineVertices[i].y == outlineVertices[j].y)
            {
                if (abs(point.y - outlineVertices[i].y) < 1e-6 &&
                    point.x >= std::min(outlineVertices[i].x, outlineVertices[j].x) &&
                    point.x <= std::max(outlineVertices[i].x, outlineVertices[j].x))
                {
                    return true;
                }
            }
            else
            {
                xIntersection = static_cast<float>((outlineVertices[j].x - outlineVertices[i].x) * (point.y - outlineVertices[i].y) / (outlineVertices[j].y - outlineVertices[i].y) + outlineVertices[i].x);
                if (abs(point.x - xIntersection) < epsilon)
                {
                    if (point.x >= std::min(outlineVertices[i].x, outlineVertices[j].x) && point.x <= std::max(outlineVertices[i].x, outlineVertices[j].x) &&
                        point.y >= std::min(outlineVertices[i].y, outlineVertices[j].y) && point.y <= std::max(outlineVertices[i].y, outlineVertices[j].y))
                    {
                        return true;
                    }
                }
            }
        }

        return false;
    }

    /**
     * Evaluates whether a given point is inside any of the shape's cutouts.
     *
     * @param point The point (Point object) in question.
     * @return true: the point is inside a cutout | false: the point is not inside a cutout.
     */
    bool isPointInsideCutout(const Point& point)
    {
        bool isInside;

        for (size_t i = 0; i < cutouts.size(); i++) // For each std::vector<Point>
        {
            isInside = false;

            for (size_t j = 0, k = cutouts[i].size() - 1; j < cutouts[i].size(); k = j++)
            {
                if (((cutouts[i][j].y > point.y) != (cutouts[i][k].y > point.y))
                    && (point.x < (cutouts[i][k].x - cutouts[i][j].x) * (point.y - cutouts[i][j].y) / (cutouts[i][k].y - cutouts[i][j].y) + cutouts[i][j].x))
                {
                    isInside = !isInside; // Even amount of intersections = OUTSIDE, Odd amount of intersections = INSIDE
                }
            }
            if (isInside) return true;
        }

        return false;
    }

    /**
     * Evaluates whether a given point is on an edge of any of the shape's cutouts.
     *
     * @param point The point (Point object) in question.
     * @param epsilon Tolerance threshold (to account for numerical precision issues). Defaults to 1e-6.
     * @return true: the point is on the edge | false: the point is not on the edge.
     */
    bool isPointOnCutoutEdge(const Point& point, float epsilon = 1e-6)
    {
        size_t numVertices = outlineVertices.size();
        float xIntersection;

        for (size_t i = 0; i < cutouts.size(); i++) // Loop through each cutout
        {
            size_t cutoutSize = cutouts[i].size();
            for (size_t j = 0, k = cutoutSize - 1; j < cutoutSize; k = j++) // Loop through each vertex of ith cutout
            {
                // Check for horizontal edge (avoid division by zero)
                if (cutouts[i][j].y == cutouts[i][k].y)
                {
                    if (abs(point.y - cutouts[i][j].y) < epsilon &&
                        point.x >= std::min(cutouts[i][j].x, cutouts[i][k].x) &&
                        point.x <= std::max(cutouts[i][j].x, cutouts[i][k].x))
                    {
                        return true;
                    }
                }
                else // If edge isn't horizontal
                {
                    // Conserve precision for comparison with tolerance threshold
                    // Perhaps cast individual operations that may lead to precision loss (divisions)?
                    xIntersection = static_cast<float>((cutouts[i][k].x - cutouts[i][j].x) * (point.y - cutouts[i][j].y) / (cutouts[i][k].y - cutouts[i][j].y) + cutouts[i][j].x);
                    if (abs(point.x - xIntersection) < 1e-6)
                    {
                        if (point.x >= std::min(cutouts[i][j].x, cutouts[i][k].x) && point.x <= std::max(cutouts[i][j].x, cutouts[i][k].x) &&
                            point.y >= std::min(cutouts[i][j].y, cutouts[i][k].y) && point.y <= std::max(cutouts[i][j].y, cutouts[i][k].y))
                        {
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }
};

int main(int argc, char* argv[])
{
    // Expect executable, text file, point x, point y
    // Throw error if arguments aren't given
    if (argc < 4) {
        std::cerr << "Please give expected parameters (text file name, point x, point y)" << std::endl;
        throw std::runtime_error("Expected text file name, point x, point y parameters.");
    }

    // Get text file name and point coordinates from command line arguments
    char* fileName = argv[1];
    int x = std::stoi(argv[2]), y = std::stoi(argv[3]);

    // Construct point from coordinates arguments
    Point point(x, y);

    // Construct shape from vertices in text file
    Shape shape = Shape::readShapeFromFile(fileName);

    // Check if the point is inside the shape
    bool inside = shape.isPointInside(point);
    if (inside)
    {
        std::cout << "Point " << point.x << ", " << point.y << " is INSIDE the shape" << std::endl;
        return 1;
    }
    else
    {
        std::cout << "Point " << point.x << ", " << point.y << " is OUTSIDE the shape" << std::endl;
        return 0;
    }
}
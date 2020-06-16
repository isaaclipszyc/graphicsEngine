#define SDL_MAIN_HANDLED
#include <ModelTriangle.h>
#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <fstream>
#include <vector>
#include <map>
#include <RayTriangleIntersection.h>

using namespace std;
using namespace glm;

#define WIDTH 640
#define HEIGHT 480

const int imageWidth = WIDTH;
const int imageHeight = HEIGHT;

int fileCount = 0;

bool saveImages = false;

vec3 cameraPos = vec3(0.f, 3.f, 2.f);
vec3 cameraDir = vec3(0.f, 0.f, 1.f);
vec3 cameraUp = vec3(0.f, 1.f, 0.f);

vec3 oldPosition = vec3(0.f, 0.f, 0.f);
mat3 oldOrientation = mat3(1.f);

mat3 cameraOrientation = mat3(1.f);

vec3 proximityLight = vec3(-0.25, 5.0, -3.0);

string currentState = "raytrace";

const float focalLength = 200.0;

vector<vec3> positions;
vector<string> colours;
vector<vector<int>> placement;
vector<vec3> colourNumber(10, vec3(5,5,5));

vector<vector<uint32_t>> texture;
vector<TexturePoint> texturePoints;

bool objAndMtlRead = false;
vector<ModelTriangle> modelTriangles;

float depthBuffer[imageWidth][imageHeight] = {};

void savePPM();
void drawLineTexture(CanvasPoint from, CanvasPoint to, TexturePoint start, TexturePoint end);
void fillTriangleTexture(CanvasTriangle points, int currentTriangle);
vec3 getNormalInverted(RayTriangleIntersection triangle);
void lookAt();
void horizontalLookAt();
void verticalLookAt();
bool nearAndFarPlaneClipping(ModelTriangle triangle);
bool backFaceCulling(ModelTriangle triangle, vec3 rayDirection);
vec3 getNormal(ModelTriangle triangle);
void initializeDepthbuffer();
void drawLine(CanvasPoint from, CanvasPoint to, uint32_t colour);
float softShadow(RayTriangleIntersection triangle, int currentTriangle);
vec3 phongShading(ModelTriangle triangle, float u, float v);
void getVertexNormals();
float gouraudShading(RayTriangleIntersection triangle, float u, float v);
void getModelTriangles();
bool inShadow(RayTriangleIntersection intersectionTriangle, int currentTriangle);
float lighting(RayTriangleIntersection triangle);
float specularLighting(RayTriangleIntersection triangle);
vec3 getNormal(RayTriangleIntersection triangle);
float proximity(vec3 point);
float angleOfIncidence(RayTriangleIntersection triangle);
void wireframe();
void rayTracing();
void fillTriangle(CanvasTriangle triangle);
vector<CanvasTriangle> orderImagePlaneTriangles(vector<CanvasTriangle> triangles);
void readObj(string filename);
void readMtl(string filename);
void fetchTriangle();
bool checkIfInside(float x, float y);
vector<vector<uint32_t>> loadImage(string filename);
void update();
void handleEvent(SDL_Event event);
vector<float> interpolate(float from, float to, float numberOfValues);
vector<vec3> interpolate3d(vec3 from, vec3 to, float numberOfValues);
vector<CanvasPoint> interpolate2d(CanvasPoint from, CanvasPoint to, float numberOfValues);
vector<TexturePoint> interpolate2d(TexturePoint from, TexturePoint to, float numberOfValues);
vector<CanvasPoint> order(CanvasTriangle triangle);


DrawingWindow window = DrawingWindow(imageWidth, imageHeight, false);

int main(int argc, char* argv[])
{
  SDL_Event event;
  if (objAndMtlRead == false) {
    readObj("cornell-box.obj");
    readMtl("cornell-box.mtl");
    getModelTriangles();
    getVertexNormals();
    objAndMtlRead = true;
  }
  
  char str[] = "texture.ppm";
  texture = loadImage(str);

  while(true)
  {
    // We MUST poll for events - otherwise the window will freeze !
    if(window.pollForInputEvents(&event)) handleEvent(event);
    update();
    window.clearPixels();

    initializeDepthbuffer();
    
    if(currentState == "wireframe") wireframe();
    if(currentState == "rasterize") fetchTriangle();
    if(currentState == "raytrace") rayTracing();

    if(saveImages) savePPM();

    // Need to render the frame at the end, or nothing actually gets shown on the screen !
    window.renderFrame();
  }
}

void initializeDepthbuffer(){
  for(int y = 0; y < imageHeight; y++){
    for(int x = 0; x < imageWidth; x++){
      depthBuffer[x][y] = std::numeric_limits<float>::infinity();
    }
  }
}

void readObj(string filename) {
  ifstream file;
  file.open(filename);

  int counter = -1;

  if(file.is_open()) {
    string line;

    while (getline(file,line)) {
      string text;
      file >> text;

      if (text == "o") {
        counter++;
      }

      if (text == "usemtl") {
        string colour;
        file >> colour;

        colours.insert(colours.end(), colour);
      }

      if (text == "v") {
        float coord1;
        float coord2;
        float coord3;

        file >> coord1;
        file >> coord2;
        file >> coord3;

        vec3 v = vec3(coord1, coord2, coord3);
        positions.push_back(v);
      }

      if (text == "vt") {
        float x;
        float y;

        file >> x;
        file >> y;

        TexturePoint t = TexturePoint(x,y);
        texturePoints.push_back(t);
      }

      if (text == "f") {
        string vertex1;
        string vertex2;
        string vertex3;
        vector<int> row;

        file >> vertex1;
        int vertex1Adjusted = stoi (vertex1);
        vertex1Adjusted = vertex1Adjusted - 1;
        row.insert(row.end(), vertex1Adjusted);

        file >> vertex2;
        int vertex2Adjusted = stoi (vertex2);
        vertex2Adjusted = vertex2Adjusted - 1;
        row.insert(row.end(), vertex2Adjusted);
        
        file >> vertex3;
        int vertex3Adjusted = stoi (vertex3);
        vertex3Adjusted = vertex3Adjusted - 1;
        row.insert(row.end(), vertex3Adjusted);

        row.insert(row.end(), counter);
        placement.push_back(row);
      }
    }
  }
}

void readMtl(string filename) {
  ifstream file;
  file.open(filename);
  string currentColour;

  if(file.is_open()) {
    string line;

    while (getline(file,line)) {

      string text;
      file >> text;

      if (text == "newmtl") {
        file >> currentColour;
      }

      if (text == "Kd") {
        float temp;
        float temp1;
        float temp2;

        file >> temp;
        file >> temp1;
        file >> temp2;

        for(int x = 0; x < (int) colours.size(); x++) {
          if(currentColour == colours[x]) {
            vec3 v = vec3(temp, temp1, temp2);
            colourNumber[x] = v;
          }
        }
      }
    }
  }
}

void fetchTriangle() {

  for(int i = 0; i < (int) modelTriangles.size(); i++) {

    ModelTriangle triangle = modelTriangles[i];

    bool clipped = nearAndFarPlaneClipping(triangle);

    if(clipped) continue;

    vector<CanvasPoint> points;

    for(int j = 0; j < 3; j++){
      vec3 vertices = cameraOrientation * (triangle.vertices[j]- cameraPos);
      CanvasPoint pixel;
      pixel.x = -(focalLength * (vertices.x / vertices.z)) + imageWidth / 2;
      pixel.y = (focalLength * (vertices.y / vertices.z)) + imageHeight / 2;
      pixel.depth =  1 / vertices.z;

      CanvasPoint canvasPoint;
      canvasPoint.x = pixel.x;
      canvasPoint.y = pixel.y;
      canvasPoint.depth = pixel.depth;

      points.push_back(canvasPoint);
    }

    vec3 colour = colourNumber[placement[i][3]];

    Colour triangleColour = Colour(round(colour[0]*255),round(colour[1]*255), round(colour[2]*255));

    // triangle on image plane
    CanvasTriangle imagePlaneTriangle = CanvasTriangle(points[0], points[1], points[2], triangleColour);
    
    fillTriangle(imagePlaneTriangle);
  }
}

bool checkIfInside(float x, float y) {
  if (round(x) > imageWidth - 1 || round(x) < 0 || round(y) > imageHeight - 1 || round(y) < 0) {
    return false;
  } else {
    return true;
  }
}

void drawLine(CanvasPoint from, CanvasPoint to, uint32_t colour) {

  float xDiff = to.x - from.x;
  float yDiff = to.y - from.y;
  float zDiff = to.depth - from.depth;

  float numberOfSteps = std::max(abs(xDiff), abs(yDiff));

  float xStep = xDiff/numberOfSteps;
  float yStep = yDiff/numberOfSteps;
  float zStep = zDiff/numberOfSteps;

  for (float i = 0.0; i < numberOfSteps; i++) {
    
    float x = from.x + (xStep * i);
    float y = from.y + (yStep * i);
    float depth = from.depth + (zStep * i);

    if(depth >= 0.0){
      continue;
    } else if(checkIfInside(x,y)){
     if (depth < depthBuffer[(int) x][(int) y]) {
        window.setPixelColour( x, y, colour);
        depthBuffer[(int) x][(int) y] = depth;
      }
    }
  }
}

void drawLineTexture(CanvasPoint from, CanvasPoint to, TexturePoint start, TexturePoint end) {

  float xDiff = to.x - from.x;
  float yDiff = to.y - from.y;
  float zDiff = to.depth - from.depth;

  float xTextDiff = end.x - start.x;
  float yTextDiff = end.y - start.y;

  float numberOfSteps = std::max(abs(xDiff), abs(yDiff));

  float xStep = xDiff / numberOfSteps;
  float yStep = yDiff / numberOfSteps;
  float zStep = zDiff / numberOfSteps;

  float xTStep = xTextDiff / numberOfSteps;
  float yTStep = yTextDiff / numberOfSteps;

  for (float i = 0.0; i < numberOfSteps; i++) {
    
    float x = from.x + (xStep * i);
    float y = from.y + (yStep * i);
    float depth = from.depth + (zStep * i);

    float xT = (start.x + (xTStep * i)) * 480;
    float yT = (start.y + (yTStep * i)) * 395;


    if(depth >= 0.0){
      continue;
    } else if(checkIfInside(x,y)){
     if (depth < depthBuffer[(int) y][(int) x]) {
        if(xT <= 480.0 && xT >= 0.0 && yT <= 395.0 && yT >= 0.0){
          cout << xT << ", " << yT << endl;
          uint32_t colour = texture[(int) round(xT)][(int) round(yT)];
          window.setPixelColour(x, y, colour);
          depthBuffer[(int) x][(int) y] = depth;
        }
      }
    }
  }
}

void fillTriangle(CanvasTriangle points) {

  vector<CanvasPoint> ordered = order(points);

  CanvasPoint topPoint = ordered[0];
  CanvasPoint midPoint = ordered[1];
  CanvasPoint bottomPoint = ordered[2];

  float red = points.colour.red;
  float green = points.colour.green;
  float blue = points.colour.blue;
  uint32_t colour = (255<<24) + (int(red)<<16) + (int(green)<<8) + int(blue);

  float prop = (midPoint.y - topPoint.y) / (bottomPoint.y - topPoint.y);
  float middlepointDepth = (prop * (bottomPoint.depth - topPoint.depth)) + topPoint.depth;

  CanvasPoint middlepoint = CanvasPoint(topPoint.x + (((midPoint.y - topPoint.y) / (bottomPoint.y - topPoint.y)) * (bottomPoint.x - topPoint.x)), midPoint.y, middlepointDepth);

  vec3 vertex1 = vec3(topPoint.x, topPoint.y, float(topPoint.depth));
  vec3 vertex2 = vec3(midPoint.x, midPoint.y, float(midPoint.depth));
  vec3 vertex3 = vec3(middlepoint.x, middlepoint.y, float(middlepoint.depth));

  drawLine(topPoint, midPoint, colour);
  drawLine(midPoint, bottomPoint, colour);
  drawLine(topPoint, bottomPoint, colour);

  //top triangle
  float diff1 = midPoint.y - topPoint.y;
  float diff2 = middlepoint.y - topPoint.y;

  vector<vec3> start = interpolate3d(vertex2, vertex1, diff1);
  vector<vec3> end = interpolate3d(vertex3, vertex1, diff2);

  for (int i = 0; i < (int) start.size(); i++) {

    CanvasPoint startPoint = CanvasPoint(start[i].x, start[i].y, start[i].z);
    CanvasPoint endPoint = CanvasPoint(end[i].x, end[i].y, end[i].z);

    drawLine(startPoint, endPoint, colour);
  }

  //bottom triangle

  vec3 vertex4 = vec3(bottomPoint.x, bottomPoint.y, float(bottomPoint.depth));

  float diff3 = bottomPoint.y - midPoint.y;
  float diff4 = bottomPoint.y - middlepoint.y;

  vector<vec3> btmStart = interpolate3d(vertex4, vertex2, diff3);
  vector<vec3> btmEnd = interpolate3d(vertex4, vertex3, diff4);

  for (int i = 0; i < (int) btmStart.size(); i++) {
    CanvasPoint startPoint = CanvasPoint(btmStart[i].x, btmStart[i].y, btmStart[i].z);
    CanvasPoint endPoint = CanvasPoint(btmEnd[i].x, btmEnd[i].y, btmEnd[i].z);

    drawLine(startPoint, endPoint, colour);
  }

  drawLine(topPoint, midPoint, colour);
  drawLine(midPoint, bottomPoint, colour);
  drawLine(topPoint, bottomPoint, colour);

}

void fillTriangleTexture(CanvasTriangle points, int currentTriangle) {

  vector<CanvasPoint> ordered = order(points);

  CanvasPoint topPoint = ordered[0];
  CanvasPoint midPoint = ordered[1];
  CanvasPoint bottomPoint = ordered[2];

  cout << texture.size() << endl;

  topPoint.texturePoint = TexturePoint(0, 0);
  midPoint.texturePoint = TexturePoint(0, 1);
  bottomPoint.texturePoint = TexturePoint(1, 1);

  float red = points.colour.red;
  float green = points.colour.green;
  float blue = points.colour.blue;
  uint32_t colour = (255<<24) + (int(red)<<16) + (int(green)<<8) + int(blue);

  float prop = (midPoint.y - topPoint.y) / (bottomPoint.y - topPoint.y);
  float middlepointDepth = (prop * (bottomPoint.depth - topPoint.depth)) + topPoint.depth;

  CanvasPoint middlepoint = CanvasPoint(topPoint.x + (((midPoint.y - topPoint.y) / (bottomPoint.y - topPoint.y)) * (bottomPoint.x - topPoint.x)), midPoint.y, middlepointDepth);

  float middlePointXText = topPoint.texturePoint.x + (midPoint.texturePoint.y - topPoint.texturePoint.y) / (bottomPoint.texturePoint.y - topPoint.texturePoint.y) * (bottomPoint.texturePoint.x - topPoint.texturePoint.x);
  float middlePointYText = midPoint.texturePoint.y;
  middlepoint.texturePoint = TexturePoint(middlePointXText, middlePointYText);
  
  vec3 vertex1 = vec3(topPoint.x, topPoint.y, float(topPoint.depth));
  vec3 vertex2 = vec3(midPoint.x, midPoint.y, float(midPoint.depth));
  vec3 vertex3 = vec3(middlepoint.x, middlepoint.y, float(middlepoint.depth));

  drawLine(topPoint, midPoint, colour);
  drawLine(midPoint, bottomPoint, colour);
  drawLine(topPoint, bottomPoint, colour);

  //top triangle
  float diff1 = midPoint.y - topPoint.y;
  float diff2 = middlepoint.y - topPoint.y;

  vector<vec3> start = interpolate3d(vertex2, vertex1, diff1);
  vector<vec3> end = interpolate3d(vertex3, vertex1, diff2);

  float textDiff1 = (midPoint.texturePoint.y - topPoint.texturePoint.y) * 395;
  float textDiff2 = (middlepoint.texturePoint.y - topPoint.texturePoint.y) * 395;

  vector<TexturePoint> textureStart = interpolate2d(midPoint.texturePoint, topPoint.texturePoint, textDiff1);
  vector<TexturePoint> textureEnd = interpolate2d(middlepoint.texturePoint, topPoint.texturePoint, textDiff2);

  for (int i = 0; i < (int) start.size(); i++) {

    CanvasPoint startPoint = CanvasPoint(start[i].x, start[i].y, start[i].z);
    CanvasPoint endPoint = CanvasPoint(end[i].x, end[i].y, end[i].z);

    TexturePoint textureStartPoint = TexturePoint(textureStart[i].x, textureStart[i].y);
    TexturePoint textureEndPoint = TexturePoint(textureEnd[i].x, textureEnd[i].y);
    drawLineTexture(startPoint, endPoint, textureStartPoint, textureEndPoint);
  }

  //bottom triangle

  vec3 vertex4 = vec3(bottomPoint.x, bottomPoint.y, float(bottomPoint.depth));

  float diff3 = bottomPoint.y - midPoint.y;
  float diff4 = bottomPoint.y - middlepoint.y;

  vector<vec3> btmStart = interpolate3d(vertex4, vertex2, diff3);
  vector<vec3> btmEnd = interpolate3d(vertex4, vertex3, diff4);

  float textDiff3 = (bottomPoint.texturePoint.y - midPoint.texturePoint.y) * 395;
  float textDiff4 = (bottomPoint.texturePoint.y - middlepoint.texturePoint.y) * 395;

  vector<TexturePoint> btmtextureStart = interpolate2d(bottomPoint.texturePoint, midPoint.texturePoint, textDiff3);
  vector<TexturePoint> btmtextureEnd = interpolate2d(bottomPoint.texturePoint, middlepoint.texturePoint, textDiff4);


  for (int i = 0; i < (int) btmStart.size(); i++) {
    CanvasPoint startPoint = CanvasPoint(btmStart[i].x, btmStart[i].y, btmStart[i].z);
    CanvasPoint endPoint = CanvasPoint(btmEnd[i].x, btmEnd[i].y, btmEnd[i].z);

    TexturePoint textureStartPoint = TexturePoint(btmtextureStart[i].x, btmtextureStart[i].y);
    TexturePoint textureEndPoint = TexturePoint(btmtextureEnd[i].x, btmtextureEnd[i].y);

    drawLineTexture(startPoint, endPoint, textureStartPoint, textureEndPoint);
  }

  drawLine(topPoint, midPoint, colour);
  drawLine(midPoint, bottomPoint, colour);
  drawLine(topPoint, bottomPoint, colour);

}

vector<CanvasPoint> order(CanvasTriangle triangle){

    CanvasPoint one = triangle.vertices[0];
    CanvasPoint two = triangle.vertices[1];
    CanvasPoint three = triangle.vertices[2];

    one.texturePoint = triangle.vertices[0].texturePoint;
    two.texturePoint = triangle.vertices[1].texturePoint;
    three.texturePoint = triangle.vertices[2].texturePoint;

    if (one.y > two.y) {
      swap(one, two);
    }
    if (two.y > three.y) {
      swap(two, three);
    }
    if (one.y > two.y) {
      swap(one, two);
    }

    vector<CanvasPoint> output;

    output.push_back(one);
    output.push_back(two);
    output.push_back(three);

    return output;
} 

void wireframe(){

  vector<CanvasTriangle> imagePlaneTriangles;

  // print face
  for(int i = 0; i < (int) modelTriangles.size(); i++) {
  
    ModelTriangle triangle = modelTriangles[i];

    bool clipped = nearAndFarPlaneClipping(triangle);

    if(clipped) continue;

    vector<CanvasPoint> points;

    for(int j = 0; j < 3; j++){

      vec3 vertices = cameraOrientation * (triangle.vertices[j] - cameraPos);
      CanvasPoint pixel;
      pixel.x = -(focalLength * (vertices.x / vertices.z)) + imageWidth / 2;
      pixel.y = (focalLength * (vertices.y / vertices.z)) + imageHeight / 2;
      pixel.depth =  1 / vertices.z;

      CanvasPoint canvasPoint;
      canvasPoint.x = pixel.x;
      canvasPoint.y = pixel.y;
      canvasPoint.depth = pixel.depth;

      points.push_back(canvasPoint);

    }
      vec3 colour = colourNumber[placement[i][3]];

      Colour triangleColour = Colour(round(colour[0]*255),round(colour[1]*255), round(colour[2]*255));

      // triangle on image plane
      CanvasTriangle imagePlaneTriangle = CanvasTriangle(points[0], points[1], points[2], triangleColour);
      imagePlaneTriangles.push_back(imagePlaneTriangle);

  }

  for(int i = 0; i < (int) imagePlaneTriangles.size(); i++){
    vector<CanvasPoint> ordered = order(imagePlaneTriangles[i]);

    CanvasPoint one = ordered[0];
    CanvasPoint two = ordered[1];
    CanvasPoint three = ordered[2];

    float red = imagePlaneTriangles[i].colour.red;
    float green = imagePlaneTriangles[i].colour.green;
    float blue = imagePlaneTriangles[i].colour.blue;
    uint32_t colour = (255<<24) + (int(red)<<16) + (int(green)<<8) + int(blue);

    drawLine(one, two, colour);
    drawLine(two, three, colour);
    drawLine(one, three, colour);
  }
}

void getModelTriangles(){
  vector<ModelTriangle> triangles;

  for(int i = 0; i < (int) placement.size(); i++) {

    vec3 p1 = positions[placement[i][0]];
    vec3 p2 = positions[placement[i][1]];
    vec3 p3 = positions[placement[i][2]]; 
    vec3 colour = colourNumber[placement[i][3]];

    Colour triangleColour = Colour(round(colour[0]*255),round(colour[1]*255), round(colour[2]*255));
    if(colour[0] == 1.0 && colour[1] == 1.0 && colour[2] == 0.0) triangleColour.name = "yellow";
    if(colour[0] == 1.0 && colour[1] == 0.0 && colour[2] == 1.0) triangleColour.name = "magenta";
    if(colour[0] == 0.700 && colour[1] == 0.700 && colour[2] == 0.700) triangleColour.name = "grey";
    if(colour[0] == 0.0 && colour[1] == 0.0 && colour[2] == 1.0) triangleColour.name = "blue";
    if(colour[0] == 0.0 && colour[1] == 1.0 && colour[2] == 0.0) {
      triangleColour.name = "green";
    }

    // triangle in 3d space, to be used for ray tracing n shit
    ModelTriangle triangle = ModelTriangle(p1, p2, p3, triangleColour);
    triangles.push_back(triangle);
  }

  modelTriangles = triangles;
}

void rayTracing(){

  RayTriangleIntersection rayIntersection;

  for(int x = 0; x < imageWidth; x++){
    for(int y = 0; y < imageHeight; y++){
      vec3 rayDirection = (cameraPos - vec3((float) ((imageWidth/2) - x), (float) (y - (imageHeight/2)), (float) focalLength)) * (cameraOrientation);

      float smallestT = numeric_limits<float>::infinity();

      float rayU = 0.0;
      float rayV = 0.0;

      int currentTriangle = 0;

      rayIntersection.intersectedTriangle.colour = Colour(0,0,0);

      for(int i = 0; i < ((int) modelTriangles.size()); i++){
        ModelTriangle triangle = modelTriangles[i];

        bool culled = backFaceCulling(triangle, rayDirection);
        bool clipped = nearAndFarPlaneClipping(triangle);

        if(culled) continue;
        if(clipped) continue;

        vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
        vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
        vec3 SPVector = cameraPos - triangle.vertices[0];
        mat3 DEMatrix(-rayDirection, e0, e1);
        vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

        float t = possibleSolution[0];
        float u = possibleSolution[1];
        float v = possibleSolution[2];

        if(t < smallestT && t >= 0.0){
          if(u >= 0.0 && u <= 1.0 && v >= 0.0 && v <= 1.0 && (u + v <= 1.0)){
            vec3 intersectPoint = cameraPos + (rayDirection * t);
            rayIntersection.intersectionPoint = intersectPoint;
            rayIntersection.intersectedTriangle = triangle;
            rayIntersection.normal = getNormal(triangle);
            smallestT = t;
            rayU = u;
            rayV = v;
            currentTriangle = i;
          }
        }
      }

      float smallestTRef = numeric_limits<float>::infinity();

      // make the blue box a mirror, well try to at least... I have spent hours trying to get this to work.
      if(rayIntersection.intersectedTriangle.colour.name == "blue"){
        vec3 normalN = rayIntersection.normal;
        vec3 reflection = (rayDirection - (2 * dot(rayDirection, normalN) * normalN));

        for(int j = 0; j < (int) modelTriangles.size(); j++){
          ModelTriangle triangleR = modelTriangles[j];
          vec3 e0N = triangleR.vertices[1] - triangleR.vertices[0];
          vec3 e1N = triangleR.vertices[2] - triangleR.vertices[0];
          vec3 SPVectorN = (rayIntersection.intersectionPoint) - triangleR.vertices[0];
          mat3 DEMatrixN(-reflection, e0N, e1N);
          vec3 possibleSolutionN = glm::inverse(DEMatrixN) * SPVectorN;

          float t = possibleSolutionN[0];
          float u = possibleSolutionN[1];
          float v = possibleSolutionN[2];

          if(t < smallestTRef && t >= 0.0){
            if(u >= 0.0 && u <= 1.0 && v >= 0.0 && v <= 1.0 && (u + v <= 1.0)){
              vec3 newintersectPoint = ((rayIntersection.intersectionPoint)) + (reflection * t);
              rayIntersection.intersectedTriangle = triangleR;
              rayIntersection.intersectionPoint = newintersectPoint;
              smallestTRef = t;
              rayU = u;
              rayV = v;
            }
          }
        }
      }

      //gouraud shading - to use comment out the next 3 lines under!
      //float gbrightness = gouraudShading(rayIntersection, rayU, rayV);
      
      //lighting using phong
      vec3 phongNormal = phongShading(rayIntersection.intersectedTriangle, rayU, rayV);
      rayIntersection.normal = phongNormal;
      float gbrightness = lighting(rayIntersection);

      bool withinShadow = inShadow(rayIntersection, currentTriangle);

      if(withinShadow == true){
        float shadow = softShadow(rayIntersection, currentTriangle);
        gbrightness =  0.1 * (1 - shadow);
      } else if(withinShadow == false && gbrightness < 0.3){
        gbrightness = 0.3;
      }

      if(gbrightness > 1.0){
        gbrightness = 1.0;
      } else if (gbrightness < 0.0){
        gbrightness = 0.0;
      }

      Colour triangleColour = rayIntersection.intersectedTriangle.colour;
      if(triangleColour.name == "blue"){
        triangleColour = Colour(0,0,0);
      }

      float red = (triangleColour.red * gbrightness);
      float green = (triangleColour.green * gbrightness);
      float blue = (triangleColour.blue * gbrightness);

      if(red > 255) red = 255;
      if(green > 255) green = 255;
      if(blue > 255) blue = 255;
      uint32_t colour = (255<<24) + ( (int) red<<16) + ((int) green<<8) + (int) blue;
      
      if (rayIntersection.distanceFromCamera < INFINITY) {
          window.setPixelColour(x, y, colour);
      }
    }
  }
}

bool nearAndFarPlaneClipping(ModelTriangle triangle){

  vec3 a = triangle.vertices[0];
  vec3 b = triangle.vertices[1];
  vec3 c = triangle.vertices[2];
  vec3 middlePoint = (a + b + c) / 3.f;

  float distanceToCamera = glm::distance(cameraPos, middlePoint);

  if(distanceToCamera <= 1.0 || distanceToCamera >= 30.0){
    return true;
  } else {
    return false;
  }

}

bool backFaceCulling(ModelTriangle triangle, vec3 rayDirection){

  vec3 normal = getNormal(triangle);

  float dotProduct = dot(normal, -rayDirection);

  if(dotProduct < 0.0){
    return true;
  } else {
    return false;
  }
}

float proximity(vec3 point){
    float distanceFromLight = glm::distance(proximityLight, point);

    float brightness = 250 / (2 * M_PI * pow(distanceFromLight,2));

    return brightness;
}

vec3 phongShading(ModelTriangle triangle, float u, float v){
  vec3 normal0 = triangle.normals[0];
  vec3 normal1 = triangle.normals[1];
  vec3 normal2 = triangle.normals[2];

  vec3 normal = ((1.f - (u + v)) * normal0) + (u * normal1) + (v * normal2);

  return normal;
}

float specularLighting(RayTriangleIntersection triangle){
  vec3 normal = triangle.normal;

  vec3 view = normalize(cameraPos - triangle.intersectionPoint);
  vec3 lightDir = normalize(triangle.intersectionPoint - proximityLight);

  vec3 reflection = normalize(lightDir - (2.0f * normal * dot(lightDir, normal)));

  float specular = pow(dot(view, reflection), 30.0f);

  if(specular < 0.0) specular = 0.0;

  return 0.4 * specular;
}

float lighting(RayTriangleIntersection triangle){

  float aoi = angleOfIncidence(triangle);

  float specular = specularLighting(triangle);

  return aoi + specular;
}

float angleOfIncidence(RayTriangleIntersection triangle){
  float brightnessProx = proximity(triangle.intersectionPoint);

  vec3 normal = triangle.normal;

  vec3 dirToLight = proximityLight - triangle.intersectionPoint;

  dirToLight = normalize(dirToLight);
  
  float angle = dot(normal, dirToLight);

  float aoi = 0.2 * (brightnessProx) +  0.5 * (angle);

  if(aoi < 0) aoi = 0;

  return 0.3 * aoi;
}

bool inShadow(RayTriangleIntersection intersectionTriangle, int currentTriangle){
  vec3  shadowRay = normalize(proximityLight - intersectionTriangle.intersectionPoint);
  float distance = glm::distance(proximityLight, intersectionTriangle.intersectionPoint);

  for(int i = 0; i < ((int) modelTriangles.size()); i++){
        ModelTriangle triangle = modelTriangles[i];
        vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
        vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
        vec3 SPVector = intersectionTriangle.intersectionPoint - triangle.vertices[0];
        mat3 DEMatrix(-shadowRay, e0, e1);
        vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

        float t = possibleSolution[0];
        float u = possibleSolution[1];
        float v = possibleSolution[2];

        if(t < distance && t > 0.001 && (currentTriangle != i)){
          if(u >= 0.0 && u <= 1.0 && v >= 0.0 && v <= 1.0 && (u + v <= 1.0)){
            return true;
          }
        }
      }
      return false;
}

void getVertexNormals(){

  vector<RayTriangleIntersection> rayTriangles;

  for(int i = 0; i < (int) modelTriangles.size(); i++){
    ModelTriangle triangle = modelTriangles[i];
    vec3 vertex1= triangle.vertices[0];
    vec3 vertex2= triangle.vertices[1];
    vec3 vertex3= triangle.vertices[2]; 

    vec3 e0 = vertex2 - vertex1;
    vec3 e1 = vertex3 - vertex1;

    vec3 normal = cross(e0, e1);

    normal = normalize(normal);

    RayTriangleIntersection rayTriangle;

    rayTriangle.normal = normal;
    rayTriangle.intersectedTriangle = triangle;

    rayTriangles.push_back(rayTriangle);
  }

  for(int i = 0; i < (int) modelTriangles.size(); i++){
    vec3 vertex0Normals = vec3(0.0,0.0,0.0);
    float count0 = 0.0;
    vec3 vertex1Normals = vec3(0.0,0.0,0.0);
    float count1 = 0.0;
    vec3 vertex2Normals = vec3(0.0,0.0,0.0);
    float count2 = 0.0;

    for(int j = 0; j < (int) rayTriangles.size(); j++){
      ModelTriangle model = modelTriangles[i];
      RayTriangleIntersection ray = rayTriangles[j];
      if(ray.intersectedTriangle.vertices[1] == model.vertices[0] || ray.intersectedTriangle.vertices[2] == model.vertices[0] || ray.intersectedTriangle.vertices[0] == model.vertices[0] ){
        vertex0Normals += ray.normal;
        count0 += 1.0;
      }
      if(ray.intersectedTriangle.vertices[0] == model.vertices[1] || ray.intersectedTriangle.vertices[2] == model.vertices[1] || ray.intersectedTriangle.vertices[1] == model.vertices[1]){
        vertex1Normals += ray.normal;
        count1 += 1.0;
      }
      if(ray.intersectedTriangle.vertices[0] == model.vertices[2] || ray.intersectedTriangle.vertices[1] == model.vertices[2] || ray.intersectedTriangle.vertices[2] == model.vertices[2]){
        vertex2Normals += ray.normal;
        count2 += 1.0;
      }
    }

    vec3 avgVertex0Normal = vertex0Normals / count0;
    vec3 avgVertex1Normal = vertex1Normals / count1;
    vec3 avgVertex2Normal = vertex2Normals / count2;

    modelTriangles[i].normals[0] = avgVertex0Normal;
    modelTriangles[i].normals[1] = avgVertex1Normal;
    modelTriangles[i].normals[2] = avgVertex2Normal;
  }  
}

vec3 getNormal(RayTriangleIntersection triangle){
  vec3 vertex1 = triangle.intersectedTriangle.vertices[0];
  vec3 vertex2 = triangle.intersectedTriangle.vertices[1];
  vec3 vertex3 = triangle.intersectedTriangle.vertices[2]; 

  vec3 e0 = vertex2 - vertex1;
  vec3 e1 = vertex3 - vertex1;

  vec3 normal = cross(e0, e1);

  normal = normalize(normal);
  
  return normal;
}

vec3 getNormalInverted(RayTriangleIntersection triangle){
  vec3 vertex1 = triangle.intersectedTriangle.vertices[2];
  vec3 vertex2 = triangle.intersectedTriangle.vertices[1];
  vec3 vertex3 = triangle.intersectedTriangle.vertices[0]; 

  vec3 e0 = vertex2 - vertex1;
  vec3 e1 = vertex3 - vertex1;

  vec3 normal = cross(e0, e1);
  normal = normalize(normal);
  
  return normal;
}

vec3 getNormal(ModelTriangle triangle){
  vec3 vertex1 = triangle.vertices[0];
  vec3 vertex2 = triangle.vertices[1];
  vec3 vertex3 = triangle.vertices[2]; 

  vec3 e0 = vertex2 - vertex1;
  vec3 e1 = vertex3 - vertex1;

  vec3 normal = cross(e0, e1);

  normal = normalize(normal);
  
  return normal;
}

float gouraudShading(RayTriangleIntersection triangle, float u, float v){

  float brightnessProx0 = proximity(triangle.intersectedTriangle.vertices[0]);
  float brightnessProx1 = proximity(triangle.intersectedTriangle.vertices[1]);
  float brightnessProx2 = proximity(triangle.intersectedTriangle.vertices[2]);

  vec3 dirToLight0 = normalize(proximityLight - triangle.intersectedTriangle.vertices[0]);
  vec3 dirToLight1 = normalize(proximityLight - triangle.intersectedTriangle.vertices[1]);
  vec3 dirToLight2 = normalize(proximityLight - triangle.intersectedTriangle.vertices[2]);

  float angle0 = dot(triangle.intersectedTriangle.normals[0], dirToLight0);
  float angle1 = dot(triangle.intersectedTriangle.normals[1], dirToLight1);
  float angle2 = dot(triangle.intersectedTriangle.normals[2], dirToLight2);

  float aoi0 = 0.2 * brightnessProx0 + 0.5 * angle0;
  float aoi1 = 0.2 * brightnessProx1 + 0.5 * angle1;
  float aoi2 = 0.2 * brightnessProx2 + 0.5 * angle2;

  if(aoi0 < 0) aoi0 = 0;
  if(aoi1 < 0) aoi1 = 0;
  if(aoi2 < 0) aoi2 = 0;

  vec3 view0 = normalize(cameraPos - triangle.intersectedTriangle.vertices[0]);
  vec3 view1 = normalize(cameraPos - triangle.intersectedTriangle.vertices[1]);
  vec3 view2 = normalize(cameraPos - triangle.intersectedTriangle.vertices[2]);

  vec3 lightDir0 = normalize(triangle.intersectedTriangle.vertices[0] - proximityLight);
  vec3 lightDir1 = normalize(triangle.intersectedTriangle.vertices[1] - proximityLight);
  vec3 lightDir2 = normalize(triangle.intersectedTriangle.vertices[2] - proximityLight);

  vec3 reflection0 = normalize(lightDir0 - (2.0f * triangle.intersectedTriangle.normals[0] * dot(lightDir0, triangle.intersectedTriangle.normals[0])));
  vec3 reflection1 = normalize(lightDir1 - (2.0f * triangle.intersectedTriangle.normals[1] * dot(lightDir1, triangle.intersectedTriangle.normals[1])));
  vec3 reflection2 = normalize(lightDir2 - (2.0f * triangle.intersectedTriangle.normals[2] * dot(lightDir2, triangle.intersectedTriangle.normals[2])));

  float specular0 = pow(dot(view0, reflection0), 1.0f);
  float specular1 = pow(dot(view1, reflection1), 1.0f);
  float specular2 = pow(dot(view2, reflection2), 1.0f);

  if(specular0 < 0.0) specular0 = 0.0;
  if(specular1 < 0.0) specular1 = 0.0;
  if(specular2 < 0.0) specular2 = 0.0;

  float brightness0 =  (0.4 * specular0) + (0.3 * aoi0);
  float brightness1 =  (0.4 * specular1) + (0.3 * aoi1);
  float brightness2 =  (0.4 * specular2) + (0.3 * aoi2);

  if(brightness0 > 1.0 ) brightness0 = 1.0;
  if(brightness0 < 0.0 ) brightness0 = 0.0;

  if(brightness1 > 1.0 ) brightness1 = 1.0;
  if(brightness1 < 0.0 ) brightness1 = 0.0;

  if(brightness2 > 1.0 ) brightness2 = 1.0;
  if(brightness2 < 0.0 ) brightness2 = 0.0;


  float approx = ((1.0 - (u + v)) * brightness0) + (u * brightness1) + (v * brightness2);

  if(approx < 0.0) approx = 0.0;

  return approx;
}

float softShadow(RayTriangleIntersection triangle, int currentTriangle){

  vec3 normal = getNormal(triangle);
  vec3 raised = triangle.intersectionPoint + (0.1f * normal);
  vec3 lowered = triangle.intersectionPoint - (0.1f * normal);

  RayTriangleIntersection raisedTriangle = RayTriangleIntersection(raised, triangle.distanceFromCamera, triangle.intersectedTriangle);

  RayTriangleIntersection loweredTriangle = RayTriangleIntersection(raised, triangle.distanceFromCamera, triangle.intersectedTriangle);

  bool above = inShadow(raisedTriangle, currentTriangle);
  bool below = inShadow(loweredTriangle, currentTriangle);

  if(above == true && below == true){
    return 0.1;
  } else if(above == false && below == false){
    return 0.0;
  } else {
    vec3 diff = raised - lowered;
    diff = diff / 10.f;
    for(int i = 1; i < 10; i++){
      vec3 gradient = lowered + (((float) i) * diff);
      RayTriangleIntersection gradientTri = RayTriangleIntersection(gradient, triangle.distanceFromCamera, triangle.intersectedTriangle);
      if(inShadow(gradientTri, currentTriangle)){
        return ( (float) (10 - i) / (float) 10.0);
      } 
    }
  }
}

vector<vector<uint32_t>> loadImage(string filename){
    ifstream reader;
    reader.open(filename, ifstream::binary);

    vector<vector<uint32_t>> image;

    if(reader.is_open()){
      string line;
      getline(reader, line);

      if(line == "P6"){
        getline(reader, line);

        if(line[0] == '#'){
          getline(reader, line);

          int whitespace = line.find(" ");

          string w = line.substr(0, whitespace);
          int width = stoi(w);

          string h = line.substr(whitespace, line.size());
          int height = stoi(h);

          getline(reader, line);

          string max = line;
          int maxVal = stoi(max);

          while((int) image.size() < height){
            vector<uint32_t> row;
            while((int) row.size() < width){
              uint8_t r = reader.get();
              uint8_t g = reader.get();
              uint8_t b = reader.get();

              float red =  (255 * (int) r) / maxVal;
              float green = (255 * (int) g) / maxVal;
              float blue = (255 * (int) b) / maxVal;

              uint32_t rgb = (255<<24) + ((int) red<<16) + ((int) green<<8) + (int) blue;
              row.push_back(rgb);
            }
            image.push_back(row);
          }
        }
      }
    }
    reader.close();

    return image;
}

void update()
{
  // Function for performing animation (shifting artifacts or moving the camera)
}

void handleEvent(SDL_Event event)
{
  if(event.type == SDL_KEYDOWN) {
    if(event.key.keysym.sym == SDLK_LEFT) {
      cameraPos = vec3(cameraPos[0] + 0.2, cameraPos[1], cameraPos[2]);
      cout << "LEFT" << endl;
    } 
    else if(event.key.keysym.sym == SDLK_RIGHT){
      cameraPos = vec3(cameraPos[0] - 0.2, cameraPos[1], cameraPos[2]);
      cout << "RIGHT" << endl;
    }
    else if(event.key.keysym.sym == SDLK_UP) {
      cameraPos = vec3(cameraPos[0], cameraPos[1] - 0.2, cameraPos[2]);
      cout << "UP" << endl;
    }
    else if(event.key.keysym.sym == SDLK_DOWN){
      cameraPos = vec3(cameraPos[0], cameraPos[1] + 0.2, cameraPos[2]);
      cout << "DOWN" << endl;
    }
    else if(event.key.keysym.sym == SDLK_f){
      cameraPos = vec3(cameraPos[0], cameraPos[1], cameraPos[2] - 0.2);
      cout << "FORWARD" << endl;
    }
    else if(event.key.keysym.sym == SDLK_b){
      cameraPos = vec3(cameraPos[0], cameraPos[1], cameraPos[2] + 0.2);
      cout << "BACK" << endl;
    }
    else if(event.key.keysym.sym == SDLK_w){
      mat3 newOrientation = cameraOrientation;
      vec3 xOrientation = vec3(1.f,0.f,0.f);
      vec3 yOrientation = vec3(0.f, cos((M_PI/180)), -sin((M_PI/180)));
      vec3 zOrientation = vec3(0.f, sin((M_PI/180)), cos((M_PI/180)));
      cameraOrientation = newOrientation * mat3(xOrientation, yOrientation, zOrientation);
      cout << "TILT UP" << endl;
    }
    else if(event.key.keysym.sym == SDLK_s){
      mat3 newOrientation = cameraOrientation;
      vec3 xOrientation = vec3(1.f,0.f,0.f);
      vec3 yOrientation = vec3(0.f, cos(-(M_PI/180)), -sin(-(M_PI/180)));
      vec3 zOrientation = vec3(0.f, sin(-(M_PI/180)), cos(-(M_PI/180)));
      cameraOrientation = newOrientation * mat3(xOrientation, yOrientation, zOrientation);
      cout << "TILT DOWN" << endl;
    }
    else if(event.key.keysym.sym == SDLK_a){
      mat3 newOrientation = cameraOrientation;
      vec3 xOrientation = vec3(cos((M_PI/180)), 0.f, sin((M_PI/180)));
      vec3 yOrientation = vec3(0.f,1.f,0.f);
      vec3 zOrientation = vec3(-sin((M_PI/180)), 0.f, cos((M_PI/180)));  
      cameraOrientation = newOrientation * mat3(xOrientation, yOrientation, zOrientation);
      cout << "PAN LEFT" << endl;
    }
    else if(event.key.keysym.sym == SDLK_d){
      mat3 newOrientation = cameraOrientation;
      vec3 xOrientation = vec3(cos(-(M_PI/180)), 0.f, sin(-(M_PI/180)));
      vec3 yOrientation = vec3(0.f,1.f,0.f);
      vec3 zOrientation = vec3(-sin(-(M_PI/180)), 0.f, cos(-(M_PI/180))); 
      cameraOrientation = newOrientation * mat3(xOrientation, yOrientation, zOrientation);
      cout << "PAN RIGHT" << endl;
    }
    else if(event.key.keysym.sym == SDLK_r){
      cameraOrientation = mat3(1.f);
      cameraPos = vec3(0.f, 3.f, 3.f);
      cout << "Reset" << endl;
    }
    else if(event.key.keysym.sym == SDLK_1){
      currentState = "wireframe";
      cout << "Wireframe" << endl;
    }
    else if(event.key.keysym.sym == SDLK_2){
      currentState = "rasterize";
      cout << "Rasterize" << endl;
    }
    else if(event.key.keysym.sym == SDLK_3){
      currentState = "raytrace";
      cout << "Raytrace" << endl;
    } 
    else if(event.key.keysym.sym == SDLK_h){
      double orbitAngle = 3.141;
      cameraPos = vec3(cameraPos.x*cos(orbitAngle*(M_PI/180)) + (cameraPos.z)*sin(orbitAngle*(M_PI/180)), 
                      cameraPos.y,
                      -cameraPos.x*sin(orbitAngle*(M_PI/180)) + (cameraPos.z)*cos(orbitAngle*(M_PI/180)));

      horizontalLookAt();
      cout << "Orbit" << endl;
    }
    else if(event.key.keysym.sym == SDLK_l){ 
      lookAt();
      cout << "Look at" << endl;
    }
    else if(event.key.keysym.sym == SDLK_p){
      saveImages = !saveImages;
    }
  }
  else if(event.type == SDL_MOUSEBUTTONDOWN) cout << "MOUSE CLICKED" << endl;
}

void savePPM(){
  if(oldPosition != cameraPos || oldOrientation != cameraOrientation){
    char filename[16];
      sprintf(filename, "video%d.ppm", fileCount);
      FILE *fp = fopen(filename, "wb");
      fprintf(fp, "P6\n%d %d\n255\n", WIDTH, HEIGHT);
      for (int j = 0; j < HEIGHT; ++j) {
        for (int i = 0; i < WIDTH; ++i) {
          uint32_t colour = window.getPixelColour(i, j);
          static unsigned char color[3];
          color[0] = (colour >> 16) & 0xff; // red
          color[1] = (colour >> 8) & 0xff; // green
          color[2] = colour  & 0xff; // blue
          fwrite(color, 1, 3, fp);
        }
      }
      fclose(fp);
      cout << "PPM saved: " << fileCount << endl;
      fileCount++;
      oldPosition = cameraPos;
      oldOrientation = cameraOrientation;
  } 
  
}

void lookAt() {
  mat3 identity = mat3(1.0f);

  double pan = atan(cameraPos.x / cameraPos.z) * (180 / M_PI);
  pan = cameraPos.z >= 0 ? pan : pan + 180;

  vec3 xOrientation = vec3(cos(pan * (M_PI / 180)), 0.f, sin(pan * (M_PI / 180)));
  vec3 yOrientation = vec3(0.f, 1.f, 0.f);
  vec3 zOrientation = vec3(-sin(pan * (M_PI / 180)), 0.f, cos(pan * (M_PI / 180))); 

  double tilt = -atan((cameraPos.y - 3) / sqrt(pow(cameraPos.x, 2) + pow(cameraPos.z, 2))) * (180 / M_PI);

  vec3 xOrientationTilt = vec3(1.f, 0.f, 0.f);
  vec3 yOrientationTilt = vec3(0, cos(tilt * (M_PI / 180)), -sin(tilt * (M_PI / 180)));
  vec3 zOrientationTilt = vec3(0, sin(tilt * (M_PI / 180)), cos(tilt * (M_PI / 180)));

  cameraOrientation = mat3(xOrientationTilt, yOrientationTilt, zOrientationTilt) * mat3(xOrientation, yOrientation, zOrientation) * identity;
}

void horizontalLookAt() {
  mat3 identity = mat3(1.0f);

  double pan = atan((cameraPos.x) / cameraPos.z) * (180 / M_PI);
  pan = cameraPos.z >= 0 ? pan : pan + 180;

  vec3 xOrientation = vec3(cos(pan * (M_PI / 180)), 0.f, sin(pan * (M_PI / 180)));
  vec3 yOrientation = vec3(0.f, 1.f, 0.f);
  vec3 zOrientation = vec3(-sin(pan * (M_PI / 180)), 0.f, cos(pan * (M_PI / 180))); 

  double tilt = atan(cameraPos.z / (cameraPos.y - 3))*(180/M_PI);
  tilt = tilt >= 0 ? tilt - 90 : tilt + 90;

  vec3 xOrientationTilt = vec3(1.f, 0.f, 0.f);
  vec3 yOrientationTilt = vec3(0, cos(tilt * (M_PI / 180)), -sin(tilt * (M_PI / 180)));
  vec3 zOrientationTilt = vec3(0, sin(tilt * (M_PI / 180)), cos(tilt * (M_PI / 180)));

  cameraOrientation = mat3(xOrientationTilt, yOrientationTilt, zOrientationTilt) * mat3(xOrientation, yOrientation, zOrientation) * identity;
}

vector<float> interpolate(float from, float to, float numberOfValues){

  float difference = to - from;
  float step = difference / numberOfValues;

  vector<float> outputs;

  for(float i = 0.0; i < numberOfValues; i++){
    float output = from + (i * step);
    if(output < 0.0){
      output = 0.0;
    }
    outputs.insert(outputs.begin(), output);
  }

  return outputs;
}

vector<vec3> interpolate3d(vec3 from, vec3 to, float numberOfValues){

  vector<vec3> outputs;

  numberOfValues = (float) glm::max(numberOfValues, 10.0f);

  float stepX = (to.x - from.x) / numberOfValues;
  float stepY = (to.y - from.y) / numberOfValues;
  float stepZ = (to.z - from.z) / numberOfValues;

  for (float i = 0.0; i < numberOfValues; i++){
    float newX = from.x + (i * stepX);
    float newY = from.y + (i * stepY);
    float newZ = from.z + (i * stepZ);
    
    vec3 newPoint = vec3(newX, newY, newZ);
    outputs.push_back(newPoint);
  }
  return outputs;
}

vector<CanvasPoint> interpolate2d(CanvasPoint from, CanvasPoint to, float numberOfValues){

  vector<CanvasPoint> outputs;

  vector<float> xvalues = interpolate(from.x, to.x, numberOfValues);
  vector<float> yvalues = interpolate(from.y, to.y, numberOfValues);
  vector<float> zvalues = interpolate(from.depth, to.depth, numberOfValues);

  for(float i = 0.0; i < numberOfValues; i++){
    CanvasPoint output = CanvasPoint(xvalues[i], yvalues[i], (double) zvalues[i]);
    outputs.insert(outputs.begin(), output);
  }
  return outputs;
}

vector<TexturePoint> interpolate2d(TexturePoint from, TexturePoint to, float numberOfValues){

  vector<TexturePoint> outputs;

  vector<float> xvalues = interpolate(from.x, to.x, numberOfValues);
  vector<float> yvalues = interpolate(from.y, to.y, numberOfValues);

  for(float i = 0.0; i < numberOfValues; i++){
    TexturePoint output = TexturePoint(round(xvalues[i]), round(yvalues[i]));
    outputs.insert(outputs.begin(), output);
  }
  return outputs;
}

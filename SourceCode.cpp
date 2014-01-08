#include <CImg.h>
#include <iostream>
#include <cstdio>
#include <vector>
#include <cmath>
#include <cstring>
#include <map>
#include <set>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <algorithm>

using namespace cimg_library;
using namespace std;

#define inf_32 (((unsigned int)1) << 31)
#define BLACK 0
#define WHITE 255
#define PI 3.141592

template <class T>
class preProcessing {
  //-- Color Filtering
public:
  CImg < T > colorFiltering2 ( CImg < T >& img ) {
    //http://www.cs.rit.edu/~ncs/color/t_convert.html
    /*======================*/
    unsigned int W = img._width, H = img._height;
    CImg < T > mask(W, H, 1, 1);
    unsigned int x, y, i;
    float hue;
    float r, g, b;
    float rgb_min, rgb_max, delta;
    for (x = 0; x < W; x++)
      for (y = 0; y < H; y++) {

        r = img(x, y, 0, 0) / 255.;
        g = img(x, y, 0, 1) / 255.;
        b = img(x, y, 0, 2) / 255.;
        if (r < g) {
          if (g < b) rgb_max = b, rgb_min = r;
          else       rgb_max = g, rgb_min = min(r, b);
        } else {
          if (r < b) rgb_max = b, rgb_min = g;
          else       rgb_max = r, rgb_min = min(g, b);
        }
        delta = rgb_max - rgb_min;

        if (delta == 0) continue;

        if (r == rgb_max) {
          hue  = (g - b);
          hue /= delta;
          if (hue < 0) hue += 6;
        } else if(g == rgb_max) {
          hue  = (b - r);
          hue /= delta;
          hue += 2;
        } else {//b == rgb_max
          hue  = (r - g);
          hue /= delta;
          hue += 4;
        }
        hue *= 60;
        if (hue < 0) hue += 360;

        if (hue < 180) continue;
        mask(x, y, 0) = BLACK;
      }
    return mask;
  }
    
  //-- Smoothing Spatial Filter
private:
  void pixSort ( T& a, T& b ) {
    if (a > b) swapIntegers(a, b);
  }
  void swapIntegers ( T& a, T &b ) {
    if (a > b) {
      a ^= b;
      b ^= a;
      a ^= b;
    }
  }
  T obtainMedian9 ( T p[] ) {
    //source : http://ndevilla.free.fr/median/median.pdf
    pixSort(p[1], p[2]); pixSort(p[4], p[5]); pixSort(p[7], p[8]);
    pixSort(p[0], p[1]); pixSort(p[3], p[4]); pixSort(p[6], p[7]);
    pixSort(p[1], p[2]); pixSort(p[4], p[5]); pixSort(p[7], p[8]);
    pixSort(p[0], p[3]); pixSort(p[5], p[8]); pixSort(p[4], p[7]);
    pixSort(p[3], p[6]); pixSort(p[1], p[4]); pixSort(p[2], p[5]);
    pixSort(p[4], p[7]); pixSort(p[4], p[2]); pixSort(p[6], p[4]);
    pixSort(p[4], p[2]); return(p[4]);
  }
public:
  CImg < T > smoothingSpatialFilter (
      CImg < T > img,
      CImg < T >& mask,
      unsigned int color_model_idx = 3
      ) {
    /*======================*/
    unsigned int W = img._width, H = img._height;
    CImg < T > ans(W, H, 1, color_model_idx, 0);
    unsigned int x, y, i, j;
    int movX[] = {-1, 0, 1, -1, 0, 1, -1, 0, 1};
    int movY[] = {-1, -1, -1, 0, 0, 0, 1, 1, 1};
    T arr[9];
    for (x = 1; x < W - 1; x++)
      for (y = 1; y < H - 1; y++) if (mask(x, y) == BLACK){
        for (j = 0; j < color_model_idx; j++) {
          for (i = 0; i < 9; i++) arr[i] = img(x + movX[i], y + movY[i], j);
          ans(x, y, j) = obtainMedian9(arr);
        }
      } else for (j = 0; j < color_model_idx; j++) ans(x, y, j) = WHITE;

    return ans;
  }
  CImg < T > rgb2gray ( CImg < T >& img ) {
    //source http://hilbert-space.de/?p=22
    //~http://www.tannerhelland.com/3643/
    /*======================*/
    unsigned int W = img._width, H = img._height;
    CImg < T > gray(W, H);
    unsigned int x, y;
    for (x = 0; x < W; x++)
      for (y = 0; y < H; y++)
        gray(x, y) = ((img(x, y, 0) * 77)  +
                      (img(x, y, 1) * 151) +
                      (img(x, y, 2) * 28)) >> 8;

    return gray;
  }

  //-- Niblack's thresholding
  void niblackThresholding (
      CImg < T >& gray_img,
      float K = -0.2,
      float C = 10,
      unsigned int wSz = 45
      ) {
    /*======================*/
    unsigned int W = gray_img._width, H = gray_img._height;
    //unsigned int X[W][H], X2[W][H]; //X, X²
    typedef vector<unsigned int> Row;
    typedef vector<Row> Matrix;
    Matrix X(W, Row(H) );
    Matrix X2(W, Row(H) );

    unsigned int x, y;
    //-- Tabla acumulativa X, X²
    for (x = 0; x < W; x++)
      for (y = 0; y < H; y++) {
        X[x][y] = gray_img(x, y);
        X2[x][y] = X[x][y] * X[x][y];
      }
    for (x = 1; x < W; x++) X[x][0] += X[x - 1][0], X2[x][0] += X2[x - 1][0];
    for (y = 1; y < H; y++) X[0][y] += X[0][y - 1], X2[0][y] += X2[0][y - 1];
    for (x = 1; x < W; x++)
      for (y = 1; y < H; y++) {
        X[x][y]  =  X[x - 1][y] +  X[x][y - 1] +  X[x][y] -  X[x - 1][y - 1];
        X2[x][y] = X2[x - 1][y] + X2[x][y - 1] + X2[x][y] - X2[x - 1][y - 1];
      }

#define get(A, x0, y0, x1, y1) (A[(x1)][(y1)] - A[(x1)][(y0) - 1] - A[(x0) - 1][(y1)] + A[(x0) - 1][(y0) - 1])
    unsigned int midWin = (wSz >> 1);
    float n = wSz * wSz;
    float mu, sigma;
    for (x = midWin + 1; x < W - midWin - 1; x++)
      for (y = midWin + 1; y < H - midWin - 1; y++) {
        mu    = get( X, x - midWin, y - midWin, x + midWin, y + midWin) / n;
        sigma = get(X2, x - midWin, y - midWin, x + midWin, y + midWin) / n;
        sigma = sqrt(sigma - mu * mu);

        gray_img(x, y) = (gray_img(x, y) > mu + K * sigma - C) ? WHITE : BLACK;
      }
#undef get
  }
  
  void merge (
      CImg < T >& ori_img,
      CImg < T >& mask
      ) {
    /*======================*/
    unsigned int W = mask._width, H = mask._height;
    unsigned int x, y, c;

    for (x = 0; x < W; x++)
      for (y = 0; y < H; y++) if (mask(x, y) == WHITE) {
        for (c = 0; c < 3; c++) ori_img(x, y, c) = WHITE;
      }
  }
};

#define __getR(x) (((x) >> 16) & 255)
#define __getG(x) (((x) >> 8)  & 255)
#define __getB(x) ( (x)        & 255)
#define sqr(x) ((x) * (x))

class UnionFind {
  vector < unsigned int > uf, sz;
public:
  void _init ( unsigned int n ) {
    uf = vector < unsigned int > (n);
    sz = vector < unsigned int > (n);
    for (unsigned int i = 0; i < n; i++) sz[uf[i] = i] = 1;
  }
  unsigned int _find ( unsigned int p ) {
    return (uf[p] == p ? p : uf[p] = _find(uf[p]));
  }
  void _union ( unsigned int p, unsigned int q ) {
    unsigned int i = _find(p), j = _find(q);
    if (i == j) return;
    if (i < j)  sz[uf[i] = j] += sz[i];
    else        sz[uf[j] = i] += sz[j];
  }
  unsigned int _size ( unsigned int p ) {
    return sz[p];
  }
  void _setParent ( unsigned int p, unsigned int parent ) {
    uf[p] = parent;
  }
  void _setSize ( unsigned int p, unsigned int size ) {
    sz[p] = size;
  }
};

template <class T>
class Mapa {
  vector < vector < pair < T, T > > > _mp;
  unsigned int MOD = 0;
public:
  void __init ( int MOD_bits = 20 ) {
    MOD = (1 << MOD_bits) - 1;
    _mp = vector < vector < pair < T, T > > > (MOD + 1, vector < pair < T, T > > ());
  }
  void __insert ( T key, T value ) {
    _mp[key & MOD].push_back(make_pair(key, value));
  }
  bool __insert_WithOut_Repetition ( T key, T value ) {
    unsigned int idx = (key & MOD);
    for (unsigned int i = 0; i < _mp[idx].size(); i++)
      if (_mp[idx][i].first == key) return false;
    _mp[idx].push_back(make_pair(key, value));
    return true;
  }
  T __get ( T key ) {
    unsigned int idx = (key & MOD);
    for (unsigned int i = 0; i < _mp[idx].size(); i++)
      if (_mp[idx][i].first == key) return _mp[idx][i].second;
    return inf_32;
  }
};

class listPixels {
  vector < unsigned int > list;
  int i;
public:
  listPixels () {
    i = 0;
    list.clear();
  }
  listPixels ( unsigned int n ) {
    i = 0;
    list.reserve(n);
  }
  void __add ( unsigned int x, unsigned int y ) {
    //--- 0 <= x, y <= 65535
    list.push_back( (x << 16) | y);
  }
  void __insertAt ( unsigned int idx, unsigned int x, unsigned int y ) {
    //--- 0 <= x, y <= 65535
    list[idx] = (x << 16) | y;
  }
#define tmp 65535
  void __get ( unsigned int idx, unsigned int& x, unsigned int& y ) {
    x = list[idx] >> 16;
    y = list[idx] & tmp;
  }
#undef tmp
  unsigned int __len () {
    return list.size();
  }
};

template <class T>
class Regions {
  unsigned int H, W;
  UnionFind uf;
  vector < int > img;

  public:
  Regions (
      unsigned int width,
      unsigned int height,
      CImg < T >& image
      ) {
    /*======================*/
    W = width;
    H = height;
    uf._init(W * H);
    img.reserve(W * H);
    for (unsigned int y, i = 0, x = 0; x < W; x++)
      for (y = 0; y < H; y++, i++)
        img[y * W + x] = ((unsigned int)image(x, y, 0) << 16) |
                         ((unsigned int)image(x, y, 1) << 8)  |
                         (image(x, y, 2));
  }

  void kMeans (
      unsigned int K,
      unsigned int d,
      unsigned int numC, vector < unsigned int >& Color,
      vector < unsigned int >& lstSize,
      vector < unsigned int >& mu,
      vector < unsigned int >& C
      ) {
    /*======================*/
    unsigned int _d, i, k;
    float dist, minDist;
    int red, green, blue;
    vector < vector < unsigned int > > sumCluster(K, vector < unsigned int > (3, 0));
    vector < unsigned int > lenCluster(K, 0);

    set < unsigned int > S;

    for (k = 0; k < K; k++) {
      do {
        i = rand() % numC;
      } while (S.find(i) != S.end());
      mu[k] = Color[i];
      S.insert(i);
    }
    S.clear();

    for (_d = 0; _d < d; _d++) {
      for (i = 0; i < numC; i++) {
        minDist = (1 << 25);
        red     = (int)__getR(Color[i]);
        green   = (int)__getG(Color[i]);
        blue    = (int)__getB(Color[i]);
        for (k = 0; k < K; k++) {
          dist = sqrt(sqr(red   - (int)__getR(mu[k])) +
                      sqr(green - (int)__getG(mu[k])) +
                      sqr(blue  - (int)__getB(mu[k])));
          if (dist < minDist) {
            minDist = dist;
            C[i]    = k;
          }
        }
        sumCluster[C[i]][0] += red   * lstSize[i];
        sumCluster[C[i]][1] += green * lstSize[i];;
        sumCluster[C[i]][2] += blue  * lstSize[i];;
        lenCluster[C[i]] += lstSize[i];
      }

      for (k = 0; k < K; k++) {
        if (lenCluster[k] == 0) lenCluster[k] = 1;
        mu[k] = ((sumCluster[k][0] / lenCluster[k]) << 16) |
                ((sumCluster[k][1] / lenCluster[k]) << 8)  |
                ((sumCluster[k][2] / lenCluster[k]));
        sumCluster[k][0] = sumCluster[k][1] = sumCluster[k][2] = 0;
        lenCluster[k]    = 0;
      }
    }
  }
  void kMeans_CostFunction (
      vector < unsigned int >& Color,
      vector < unsigned int >& lstSize,
      vector < unsigned int >& opt_mu,
      vector < unsigned int >& opt_C,
      unsigned int K = 3,
      unsigned int t = 5,
      unsigned int d = 100
      ) {
    /*======================*/
    srand (time(NULL));
    unsigned int _t, i, numC = Color.size();
    float J, minJ = (1 << 25);
    vector < unsigned int > mu (K);
    vector < unsigned int > C (numC);

    for (_t = 0; _t < t; _t++) {
      kMeans(K, d, numC, Color, lstSize, mu, C);
      J = 0;
      for (i = 0; i < numC; i++) {
        J += sqrt(sqr((int)__getR(Color[i]) - (int)__getR(mu[C[i]])) +
                  sqr((int)__getG(Color[i]) - (int)__getG(mu[C[i]])) +
                  sqr((int)__getB(Color[i]) - (int)__getB(mu[C[i]])));
      }
      J = J / numC;
      cout << "J : " << J  << endl;
      if (J < minJ) {
        minJ = J;
        swap(mu, opt_mu);
        swap(C, opt_C);
      }
    }
    cout << "Jmin : " << " " << minJ << endl;
  }

  void compressWithMask (
      vector < T >& params,
      CImg < T >&   image,
      CImg < T >&   mask,
      char * *      argv
      ) {
    /*======================*/
    unsigned int x, y, c, mv, i, parent, next_parent;
    bool can;
    unsigned int mX[] = {0, 1, 1}, mY[] = {1, 0, 1}, lenM = 3;
    unsigned int current, next;
    unsigned int tmp;

    listPixels validPixels, validBorderPixelsX, validBorderPixelsY;
    unsigned int lenValidPixels, lenValidBorderPixelsX, lenValidBorderPixelsY;

    for (x = 0; x < W - 1; x++)
      for (y = 0; y < H - 1; y++) 
        if(mask(x, y) == BLACK) validPixels.__add(x, y);
    for (x = 0; x < W - 1; x++)
      if (mask(x, H - 1) == BLACK) validBorderPixelsX.__add(x, H - 1);
    for (y = 0; y < H - 1; y++)
      if (mask(W - 1, y) == BLACK) validBorderPixelsY.__add(W - 1, y);

    lenValidPixels = validPixels.__len();
    lenValidBorderPixelsX = validBorderPixelsX.__len();
    lenValidBorderPixelsY = validBorderPixelsY.__len();

    for (mv = 0; mv < lenM; mv++)
      for (i = 0; i < lenValidPixels; i++) {
        validPixels.__get(i, x, y);
        if (mask(x + mX[mv], y + mY[mv]) == BLACK) {
          parent      = uf._find(y * W + x);
          next_parent = (y + mY[mv]) * W + (x + mX[mv]);
          if (parent == next_parent) continue;
          current = img[parent];
          next    = img[next_parent];
          can = abs((char)__getR(current) - (char)__getR(next)) <= params[0] &&
                abs((char)__getG(current) - (char)__getG(next)) <= params[1] &&
                abs((char)__getB(current) - (char)__getB(next)) <= params[2];

          if (can)
            uf._union(parent, next_parent);
        }
      }
    //--
    mv = 0;
    for (i = 0; i < lenValidBorderPixelsY - 1; i++) {
      validBorderPixelsY.__get(i, x, y);
      if (mask(x + mX[mv], y + mY[mv]) == BLACK) {
        parent      = uf._find(y * W + x);
        next_parent = (y + mY[mv]) * W + (x + mX[mv]);
        if (parent == next_parent) continue;
        current = img[parent];
        next    = img[next_parent];
        can = abs((char)__getR(current) - (char)__getR(next)) <= params[0] &&
              abs((char)__getG(current) - (char)__getG(next)) <= params[1] &&
              abs((char)__getB(current) - (char)__getB(next)) <= params[2];
        if (can) uf._union(parent, next_parent);
      }
    }
    validBorderPixelsY.__get(lenValidBorderPixelsY - 1, x, y);
    if (mask(x + mX[mv], y + mY[mv]) == BLACK) {
      parent      = uf._find(y * W + x);
      next_parent = (y + mY[mv]) * W + (x + mX[mv]);
      current = img[parent];
      next    = img[next_parent];
      can = abs((char)__getR(current) - (char)__getR(next)) <= params[0] &&
            abs((char)__getG(current) - (char)__getG(next)) <= params[1] &&
            abs((char)__getB(current) - (char)__getB(next)) <= params[2];
      if (can) uf._union(parent, next_parent);
    }
    mv = 1;
    for (i = 0; i < lenValidBorderPixelsX - 1; i++) {
      validBorderPixelsX.__get(i, x, y);
      if (mask(x + mX[mv], y + mY[mv]) == BLACK) {
        parent      = uf._find(y * W + x);
        next_parent = (y + mY[mv]) * W + (x + mX[mv]);
        if (parent == next_parent) continue;
        current = img[parent];
        next    = img[next_parent];
        can = abs((char)__getR(current) - (char)__getR(next)) <= params[0] &&
              abs((char)__getG(current) - (char)__getG(next)) <= params[1] &&
              abs((char)__getB(current) - (char)__getB(next)) <= params[2];
        if (can) uf._union(parent, next_parent);
      }
    }
    validBorderPixelsX.__get(lenValidBorderPixelsX - 1, x, y);
    if (mask(x + mX[mv], y + mY[mv]) == BLACK) {
      parent      = uf._find(y * W + x);
      next_parent = (y + mY[mv]) * W + (x + mX[mv]);
      current = img[parent];
      next    = img[next_parent];
      can = abs((char)__getR(current) - (char)__getR(next)) <= params[0] &&
            abs((char)__getG(current) - (char)__getG(next)) <= params[1] &&
            abs((char)__getB(current) - (char)__getB(next)) <= params[2];
      if (can) uf._union(parent, next_parent);
    }

    vector < unsigned int > lstParent, lstSize;
    Mapa < unsigned int > m_idxParent;
    m_idxParent.__init();
    unsigned int numComp = 0;

    for (i = 0; i < lenValidPixels; i++) {
      validPixels.__get(i, x, y);
      parent = uf._find(y * W + x);
      if (parent == y * W + x) {
        lstParent.push_back(parent);
        lstSize.push_back(uf._size(parent));
        m_idxParent.__insert(parent, numComp);
        numComp++;
      }
    }
    for (i = 0; i < lenValidBorderPixelsX; i++) {
      validBorderPixelsX.__get(i, x, y);
      parent = uf._find(y * W + x);
      if (parent == y * W + x) {
        lstParent.push_back(parent);
        lstSize.push_back(uf._size(parent));
        m_idxParent.__insert(parent, numComp);
        numComp++;
      }
    }
    for (i = 0; i < lenValidBorderPixelsY; i++) {
      validBorderPixelsY.__get(i, x, y);
      parent = uf._find(y * W + x);
      if (parent == y * W + x) {
        lstParent.push_back(parent);
        lstSize.push_back(uf._size(parent));
        m_idxParent.__insert(parent, numComp);
        numComp++;
      }
    }

    vector < vector < unsigned int > > tmpSum ( numComp, vector < unsigned int > (3, 0));
    for (i = 0; i < lenValidPixels; i++) {
      validPixels.__get(i, x, y);
      current = y * W + x;
      parent  = uf._find(current);
      c   = img[current];
      tmp = m_idxParent.__get(parent);
      tmpSum[tmp][0] += __getR(c);
      tmpSum[tmp][1] += __getG(c);
      tmpSum[tmp][2] += __getB(c);
    }

    vector < unsigned int > lstColorParent (numComp);
    for (unsigned int sz, i = 0; i < numComp; i++) {
      sz = lstSize[i];
      lstColorParent[i] = ((tmpSum[i][0] / sz) << 16) | 
                          ((tmpSum[i][1] / sz) << 8)  |
                          ((tmpSum[i][2] / sz));
    }

    unsigned int K = 3;
    vector < unsigned int > mu(K), C(numComp);
    kMeans_CostFunction(lstColorParent, lstSize, mu, C, K);
    for (x = 0; x < K; x++) 
      cout << "color " << (x + 1) << " " 
           << __getR(mu[x]) << " " 
           << __getG(mu[x]) << " " 
           << __getB(mu[x]) << endl;
    
    char path [100];
    string str;

    vector < pair < unsigned int, unsigned int > > tmptmp (3);
    vector < unsigned int > Color(3, WHITE);
    unsigned int grayscale_cluster1 = ((__getR(mu[0]) * 77)  + 
                                       (__getG(mu[0]) * 151) + 
                                       (__getB(mu[0]) * 28)) >> 8;
    unsigned int grayscale_cluster2 = ((__getR(mu[1]) * 77)  + 
                                       (__getG(mu[1]) * 151) +
                                       (__getB(mu[1]) * 28)) >> 8;
    unsigned int grayscale_cluster3 = ((__getR(mu[2]) * 77)  +
                                       (__getG(mu[2]) * 151) + 
                                       (__getB(mu[2]) * 28)) >> 8;
    tmptmp[0] = make_pair(grayscale_cluster1, 0);
    tmptmp[1] = make_pair(grayscale_cluster2, 1);
    tmptmp[2] = make_pair(grayscale_cluster3, 2);
    sort(tmptmp.begin(), tmptmp.end());

    Color[tmptmp[0].second] = Color[tmptmp[1].second] = BLACK;
    CImg < T > ans(W, H);
    filterComponents(mask, validPixels, lenValidPixels, validBorderPixelsX, lenValidBorderPixelsX, validBorderPixelsY, lenValidBorderPixelsY, m_idxParent, Color, C, ans);
    str = ("image_01:02_" + std::string(argv[1]));
    i = 0; for (i = 0; i < str.size(); i++) path[i] = str[i]; path[i] = '\0';
    ans.save(path);
  }
  
  void filterComponents (
      CImg < T >&   mask,
      listPixels&   validPixels,
      unsigned int& lenValidPixels,
      listPixels&   validBorderPixelsX,
      unsigned int& lenValidBorderPixelsX,
      listPixels&   validBorderPixelsY,
      unsigned int& lenValidBorderPixelsY,
      Mapa < unsigned int >&   m_idxParent,
      vector < unsigned int >& color,
      vector < unsigned int >& Cluster,
      CImg < T >&   answer
      ) {
    //--- Calcular dimensiones de los BoundingBoxes
    unsigned int x, y, xx, yy, i, j, mv, parent, next_parent, tmp;
    unsigned int mX[] = {0, 1}, mY[] = {1, 0}, lenM = 2;
    Mapa < unsigned int > m_Parent;
    vector < unsigned int > im_Parent;
    m_Parent.__init(18);
    
    for (j = i = 0; i < Cluster.size(); i++) {
      if (color[Cluster[i]] == BLACK) {
        m_Parent.__insert(i, j);
        im_Parent.push_back(i);
        j++;
      }
    }

    UnionFind uf_FilterComponents;
    uf_FilterComponents._init(j);
    for (mv = 0; mv < lenM; mv++)
      for (i = 0; i < lenValidPixels; i++) {
        validPixels.__get(i, x, y);
        parent = m_idxParent.__get(uf._find(y * W + x));
        if (color[Cluster[parent]] != BLACK) continue;
        if (mask(x + mX[mv], y + mY[mv]) != BLACK) continue;
        next_parent = m_idxParent.__get(uf._find((y + mY[mv]) * W + (x + mX[mv])));
        if (color[Cluster[next_parent]] == BLACK) {
          uf_FilterComponents._union(m_Parent.__get(parent), m_Parent.__get(next_parent));
        }
      }
    unsigned int numComp = 0;
    vector < unsigned int > m_idxParentComponents (j);
    for (i = 0; i < j; i++) if (i == uf_FilterComponents._find(i)) {
      m_idxParentComponents[i] = numComp++;
    }

    //--- Calcular los boundingBoxes
    vector < pair < pair < unsigned int, unsigned int >, pair < unsigned int, unsigned int > > > lstBoundingBoxes (numComp, make_pair(make_pair(inf_32, inf_32), make_pair(0, 0)));
    vector < vector < unsigned int > > lstCompBoundingBoxes(numComp, vector < unsigned int >());

    for (i = 0; i < lenValidPixels; i++) {
      validPixels.__get(i, x, y);
      if (color[Cluster[parent = m_idxParent.__get(uf._find(y * W + x))]] != BLACK) continue;
      parent = m_idxParentComponents[ uf_FilterComponents._find(m_Parent.__get(parent)) ];
      x *= 10; //-- para evitar operaciones con punto flotante al determinar el centro de cada boundingBox
      y *= 10;
      if (x < lstBoundingBoxes[parent].first.first)   lstBoundingBoxes[parent].first.first   = x;
      if (x > lstBoundingBoxes[parent].second.first)  lstBoundingBoxes[parent].second.first  = x;
      if (y < lstBoundingBoxes[parent].first.second)  lstBoundingBoxes[parent].first.second  = y;
      if (y > lstBoundingBoxes[parent].second.second) lstBoundingBoxes[parent].second.second = y;
      lstCompBoundingBoxes[parent].push_back(i);
    }

    //--- Calcular los centros de cada BoundingBox & validarlos
    float _width, _height, area, increment, theta;
    float xmin,xmax,ymin,ymax,xtemp,ytemp,ltemp,wtemp;
    vector < pair < unsigned int, unsigned int > > lstCenterBoundingBoxes (numComp);
    vector < bool > lstFinalValidComponents (numComp, true);
    set < unsigned int > tmpX, tmpY;
    /******/
    for (i = 0; i < numComp; i++) {
      if (lstCompBoundingBoxes[i].size() < 50) {
        lstFinalValidComponents[i] = false;
        continue;
      }
    }
    /******/
    for (i = 0; i < numComp; i++) if (lstFinalValidComponents[i]) {
      if (abs(lstBoundingBoxes[i].second.first - lstBoundingBoxes[i].first.first) < 15 || abs(lstBoundingBoxes[i].second.second - lstBoundingBoxes[i].first.second) < 15) {
        lstFinalValidComponents[i] = false;
        continue;
      }
      //-- Check the aspect ratio between width & height
      area = (lstBoundingBoxes[i].second.first - lstBoundingBoxes[i].first.first) * (lstBoundingBoxes[i].second.second - lstBoundingBoxes[i].first.second);
      increment = 1./36.;
      for (theta = increment * PI; theta < PI / 2.0; theta += increment * PI) {
        xmin = ymin = 1000000;
        xmax = ymax = 0;
        for (j = 0; j < lstCompBoundingBoxes[i].size(); j++) {
          validPixels.__get(lstCompBoundingBoxes[i][j], x, y);
          xtemp = x * cos(theta) + y * -sin(theta);
          ytemp = x * sin(theta) + y *  cos(theta);
          xmin = min(xtemp,xmin);
          xmax = max(xtemp,xmax);
          ymin = min(ytemp,ymin);
          ymax = max(ytemp,ymax);
        }
        ltemp = xmax - xmin + 1;
        wtemp = ymax - ymin + 1;
        if (ltemp * wtemp < area) {
          area = ltemp * wtemp;
          _height = ltemp;
          _width = wtemp;
        }
      }
      // check if the aspect ratio is between 1/10 and 10
      if (_height / _width < 1./3. || _height / _width > 3.) {
        lstFinalValidComponents[i] = false;
        continue;
      }

      x = (lstBoundingBoxes[i].first.first  + lstBoundingBoxes[i].second.first)  >> 1;
      y = (lstBoundingBoxes[i].first.second + lstBoundingBoxes[i].second.second) >> 1;

      lstCenterBoundingBoxes[i] = make_pair(x, y);
      tmpX.insert(x);
      tmpY.insert(y);
    }

    //--- Compresion de coordenadas, tabla aditiva
    vector < unsigned int > CoordinateX, CoordinateY;
    Mapa < unsigned int > iCoordinateX, iCoordinateY;
    iCoordinateX.__init(18);
    iCoordinateY.__init(18);
    i = 0; for (set < unsigned int > :: iterator it = tmpX.begin(); it != tmpX.end(); it++) CoordinateX.push_back(*it), iCoordinateX.__insert(*it, i++);
    i = 0; for (set < unsigned int > :: iterator it = tmpY.begin(); it != tmpY.end(); it++) CoordinateY.push_back(*it), iCoordinateY.__insert(*it, i++);
    tmpX.clear(); tmpY.clear();

    typedef vector<unsigned int> Row;
    typedef vector<Row> Matrix;
    unsigned int tmpW = CoordinateX.size();
    unsigned int tmpH = CoordinateY.size();
    Matrix T_acum(tmpW + 1, Row(tmpH + 1, 0));
    for (i = 0; i < numComp; i++) {
      if (!lstFinalValidComponents[i]) continue;
      x = iCoordinateX.__get(lstCenterBoundingBoxes[i].first);
      y = iCoordinateY.__get(lstCenterBoundingBoxes[i].second);
      T_acum[x + 1][y + 1] = 1;
    }
    for (x = 2; x <= tmpW; x++) T_acum[x][0] += T_acum[x - 1][0];
    for (y = 2; y <= tmpH; y++) T_acum[0][y] += T_acum[0][y - 1];
    for (x = 2; x <= tmpW; x++) for (y = 2; y <= tmpH; y++)
      T_acum[x][y] = T_acum[x - 1][y] + T_acum[x][y - 1] + T_acum[x][y] - T_acum[x - 1][y - 1];

#define get(A, x0, y0, x1, y1) (A[(x1)][(y1)] - A[(x1)][(y0) - 1] - A[(x0) - 1][(y1)] + A[(x0) - 1][(y0) - 1])
    //--- validar componentes
    answer = answer.fill(WHITE);
    for (i = 0; i < numComp; i++) {
      if (!lstFinalValidComponents[i]) continue;
      x  = (unsigned int)(lower_bound(CoordinateX.begin(), CoordinateX.end(), lstBoundingBoxes[i].first.first)   - CoordinateX.begin());
      y  = (unsigned int)(lower_bound(CoordinateY.begin(), CoordinateY.end(), lstBoundingBoxes[i].first.second)  - CoordinateY.begin());
      xx = (unsigned int)(upper_bound(CoordinateX.begin(), CoordinateX.end(), lstBoundingBoxes[i].second.first)  - CoordinateX.begin()) - 1;
      yy = (unsigned int)(upper_bound(CoordinateY.begin(), CoordinateY.end(), lstBoundingBoxes[i].second.second) - CoordinateY.begin()) - 1;
      //-- If the number of components inner one component is greater than 2, remove it as a component
      if (get(T_acum, x + 1, y + 1, xx + 1, yy + 1) > 2) {
        lstFinalValidComponents[i] = false;
        continue;
      }
      //-- componente valida
      for (j = 0; j < lstCompBoundingBoxes[i].size(); j++) {
        validPixels.__get(lstCompBoundingBoxes[i][j], x, y);
        answer(x, y) = BLACK;
      }
    }
#undef get
  }
};

int main(int argc, char * * argv) {
  //-- ./segmentacionKM /home/user/Pictures/FotosCalles/001.jpg
  typedef unsigned char color_type;

  char path [100];
  string str = (std::string(argv[1]));
  int i = 0; for (i = 0; i < str.size(); i++) path[i] = str[i]; path[i] = '\0';
  cout << str << endl;

  //.. Open the image file, get width and height values
  CImg < color_type > src(path);
  //CImgDisplay main_disp2(src, "Original");

  preProcessing < color_type > pre;

  CImg < color_type > mask = pre.colorFiltering2(src);
  //CImgDisplay main_disp3(mask, "Filt01");

  CImg  < color_type > gray = pre.rgb2gray(src);
  //CImgDisplay main_disp4(gray, "Filt02");

  CImg < color_type > filt02 = pre.smoothingSpatialFilter(gray, mask, 1);
  //CImgDisplay main_disp5(filt02, "Filt03");
  free(mask);
  free(gray);

  pre.niblackThresholding(filt02);
  //CImgDisplay main_disp6(filt02, "Filt04");

  pre.merge(src, filt02);
  //CImgDisplay main_disp7(filt02, "Filt05");
  

  //--Disjoint sets
  color_type P[] = {5, 5, 5};
  vector < color_type > params;
  params.push_back(P[0]);
  params.push_back(P[1]);
  params.push_back(P[2]);

  unsigned int width  = mask.width();
  unsigned int height = mask.height();
  Regions < color_type > op (width, height, src);
  op.compressWithMask(params, src, filt02, argv);
  //CImgDisplay main_disp8(src, "Result");
  return 0;
}

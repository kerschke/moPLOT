#include <Rcpp.h>
#include <cmath>
#include <iostream>
#include <math.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector compute3DcrossProductCPP(NumericVector y, NumericVector z) {
  // helper function that computes the cross product between two 3D vectors
  NumericVector x(3, 0.0);
  x(0) = y(1) * z(2) - y(2) * z(1);
  x(1) = y(2) * z(0) - y(0) * z(2);
  x(2) = y(0) * z(1) - y(1) * z(0);
  return x;
}

// [[Rcpp::export]]
double computeVectorLengthCPP(NumericVector vec) {
  // computes length of the vector
  return sqrt(sum(vec * vec)) ;
}

// [[Rcpp::export]]
NumericVector normalizeVectorCPP(NumericVector vec, double prec) {
  // normalizes vector vec by its length
  NumericVector vec2 = ifelse(abs(vec) < prec, 0, vec);
  if (is_false(all(vec2 == 0))) {
    vec2 = vec2 / computeVectorLengthCPP(vec2);
  }
  return vec2 ;
}

// [[Rcpp::export]]
NumericMatrix normalizeMatrixRowsCPP(NumericMatrix mat, double prec) {
  for (int i = 0; i < mat.nrow(); i++) {
    mat(i,_) = normalizeVectorCPP(mat(i,_), prec);
  }

  return mat;
}

// [[Rcpp::export]]
double computeAngleCPP(NumericVector vec1, NumericVector vec2, double prec) {
  // computes angle spanned by vectors vec1 and vec2
  // prec = max absolute difference (per element) between the vectors

  double angle = 0;
  if (max(abs(vec1 + vec2)) < prec) {
    // if vec1 and vec2 show in opposite directions
    angle = 180;
  } else {
    if (max(abs(vec1 - vec2)) >= prec) {
      // if vec1 and vec2 are (almost) identical
      double enumerator = sum(vec1 * vec2);
      double denominator = computeVectorLengthCPP(vec1) * computeVectorLengthCPP(vec2);
      if (enumerator < denominator) {
        // the enumerator can only be larger than the denominator, in case of numerical imprecision
        double rad = acos(enumerator / denominator);
        double pi = atan(1) * 4;
        angle = 180 * rad / pi;
      }
    }
  }
  return angle;
}

// [[Rcpp::export]]
IntegerVector findNextCellCPP(NumericVector gradient) {
  // which is the next cell one has to move to (based on the high-dimensional gradient)?

  int d = gradient.size();
  IntegerVector direction(d, 0.0);
  
  NumericVector normalized_gradient = normalizeVectorCPP(gradient, 0);

  // stretch gradient vector such that its longest component touches the boundary
  // of [-3, 3]^d; for each component with a (stretched) value >= 1 go in positive
  // direction, for all components with a value <= -1 go in the negative direction
  for (int j = 0; j < d; j++) {
    
    // Threshold corresponds to an angle of 22.5 degrees in a given direction
    double threshold = sin((22.5 / 360) * 2 * M_PI);
    
    if (normalized_gradient(j) < -threshold) {
      direction(j) = -1;
    } else if (normalized_gradient(j) > threshold) {
      direction(j) = 1;
    } else {
      direction(j) = 0;
    }
  }
  return direction;
}

// [[Rcpp::export]]
int convertIndices2CellIDCPP(IntegerVector indices, IntegerVector dims) {
  // convert index per dimension [(1, ..., rows), (1, ..., columns), ...] to cell ID (1, ..., prod(dims))
  int cellID = -1;
  if (is_true(any(indices > dims)) | is_true(any(indices < 1))) {
    // if either the row- or columnIndex are located outside the boundaries, return -1 as ID
    cellID = -1;
  } else {
    int d = dims.size();
    std::vector<int> cumcells(d);
    cumcells[0] = 1;
    for (int j = 1; j < d; j++) {
      cumcells[j] = cumcells[j - 1] * dims[j - 1];
    }
    cellID = 1;
    for (int j = 0; j < d; j++) {
      cellID += (indices[j] - 1) * cumcells[j];
    }
  }
  return cellID;
}

// [[Rcpp::export]]
IntegerVector convertCellID2IndicesCPP(int cellID, IntegerVector dims) {
  // convert cellID (1, ..., prod(dims)) to index per dimension [(1, ..., rows), (1, ..., columns), ...]
  int d = dims.size();
  IntegerVector indexVector(d, -1);

  IntegerVector cumcells(d, 1);
  cumcells[0] = dims[0];
  for (int j = 1; j < d; j++) {
    cumcells[j] = cumcells[j - 1] * dims[j];
  }

  if ((cellID <= cumcells[d - 1]) & (cellID >= 1)) {
    indexVector[0] = ((cellID - 1) % cumcells[0]) + 1;
    for (int j = 1; j < d; j++) {
      int tmp = floor((cellID - 1) / cumcells[j - 1]);
      indexVector[j] = (tmp % dims[j]) + 1;
    }
  }
  return indexVector;
}

double relativeVectorOrientation(NumericVector a, NumericVector b) {
  // > 0 --> b right of a
  // < 0 --> b left of a
  // == 0 --> (anti-) parallel
  return (a(0) * -b(1)) + (a(1) * b(0));
}

// [[Rcpp::export]]
IntegerMatrix getNeighbourhood(int d, bool include_diagonals) {
  if (include_diagonals) {
    Rcpp::Function expandGrid("expand.grid");

    IntegerVector v = IntegerVector::create(-1, 0, 1);

    DataFrame deltasDF;

    if (d == 2) {
      deltasDF = expandGrid(v, v);
    } else if (d == 3) {
      deltasDF = expandGrid(v, v, v);
    } else {
      warning("Neighbourhood could not be created, d not in {2,3}!");
    }
    return internal::convert_using_rfunction(deltasDF, "as.matrix");
  } else {
    IntegerMatrix deltas(3 * d, d);
    int cnt = 0;
    for (int d_iter = 0; d_iter < d; d_iter++) {
      for (int step = -1; step <= 1; step++) {
        IntegerVector v(d, 0);
        v(d_iter) = step;
        deltas(cnt++, _) = v;
      }
    }
    return deltas;
  }
}

// [[Rcpp::export]]
NumericMatrix imputeBoundary(NumericMatrix moGradMat, List gradMatList, IntegerVector dims) {
  int p = gradMatList.length();
  int d = dims.size();
  int n = moGradMat.nrow();

  NumericMatrix moGrad = clone(moGradMat);

  IntegerVector indices;
  NumericVector newGradient;

  std::vector<NumericMatrix> gradMat(p);

  for (int i = 0; i < gradMatList.length(); i++) {
    gradMat[i] = as<NumericMatrix>(gradMatList(i));
  }

  for (int id = 1; id <= n; id++) {
    indices = convertCellID2IndicesCPP(id, dims);

    if (is_true(any(indices == 1)) || is_true(any(indices == dims))) {
      // id is at a border

      // for each dimension, if the MOG is leaving the legal area,
      // see if setting that particular dimension to zero results in a
      // vector consistent with the single-objective gradients

      for (int d_iter = 0; d_iter < d; d_iter++) {
        if ((indices(d_iter) == 1            && moGrad(id-1,d_iter) < 0) || // lower bound and negative-going component
            (indices(d_iter) == dims(d_iter) && moGrad(id-1,d_iter) > 0))   // upper bound and positive-going component
        {
          // test if setting this dimension to zero is still consistent with the
          // single-objective gradients
          double length = computeVectorLengthCPP(moGrad(id-1,_));

          newGradient = moGrad(id-1,_);
          newGradient(d_iter) = 0;

          bool legal = true;

          for (int f_id = 0; f_id < p; f_id++) {
            if (sum(gradMat[f_id](id-1,_) * newGradient) < 0) {
              legal = false;
              break;
            }
          }

          if (legal) {
            if (is_false(all(newGradient == 0))) {
              moGrad(id-1,_) = newGradient / computeVectorLengthCPP(newGradient) * length;
            } else {
              moGrad(id-1,_) = newGradient;
            }
          }
        }
      }
    }
  }

  return moGrad;
}

int isCritical(std::vector<NumericVector> vectors) {
  int d = vectors[0].size();
  
  LogicalVector positive(d, false);
  LogicalVector negative(d, false);

  for (int id = 0; id < vectors.size(); id++) {
    positive = positive | (vectors[id] >= 0);
    negative = negative | (vectors[id] <= 0);
  }

  if (is_true(any(negative & !positive)) ||
      is_true(any(positive & !negative))) {
    // cannot be critical
    return -1;
  }
  
  if (d == 2) {
    // for each vector ...
    for (int vID = 0; vID < vectors.size(); vID++) {
      bool pos = false;
      bool neg = false;
      bool zero = false;
      
      // for each other vector ...
      for (int compID = 0; compID < vectors.size(); compID++) {
        if (is_true(all(vectors[vID] == vectors[compID]))) continue;
        
        double orientation = relativeVectorOrientation(vectors[vID], vectors[compID]);
        
        // TODO >=,<= or >, < ?!
        if (orientation > 0) {
          pos = true;
        }
        if (orientation < 0) {
          neg = true;
        }
        if (orientation == 0 && (sum(vectors[vID] * vectors[compID]) <= 0)) {
          zero = true;
        }
      }
      
      if (!(pos && neg)) {
        if (zero) {
          // border case
          return 0;
        } else {
          // not critical
          return -1;
        }
      }
    }
  } else if (d == 3) {
    // for each potential plane ...
    for (int vID = 0; vID < vectors.size(); vID++) {
      for (int vID2 = 0; vID2 < vectors.size(); vID2++) {
        if (is_true(all(vectors[vID] == vectors[vID2]))) continue;
        
        bool pos = false;
        bool neg = false;
        bool zero = false;
        
        // compute the normal vector
        NumericVector normal_vector = compute3DcrossProductCPP(vectors[vID], vectors[vID2]);
        
        // if normal_vector == 0, we cannot do reasonable checks here
        if (is_true(all(normal_vector == 0))) continue;
        
        // for each other vector ...
        for (int compID = 0; compID < vectors.size(); compID++) {
          if (is_true(all(vectors[vID] == vectors[compID]))) continue;
          if (is_true(all(vectors[vID2] == vectors[compID]))) continue;
          
          double orientation = sum(normal_vector * vectors[compID]);
          
          // TODO >=,<= or >, < ?!
          if (orientation > 0) {
            pos = true;
          }
          if (orientation < 0) {
            neg = true;
          }
          if (orientation == 0) {
            zero = true;
          }
        }
        
        if (!(pos && neg)) {
          if (zero) {
            if (!pos && !neg) {
              // degenerate critical
              return 1;
            } else {
              return 0;
            }
          } else {
            // not critical
            return -1;
          }
        }
      }
    }
    
  }

  // critical
  return 1;
}

std::vector<std::vector<int>> getCellSimplices(int dimension) {
  if (dimension == 2) {
    std::vector<std::vector<int>> simplices {
      {
        0b00,
        0b01,
        0b11
      },{
        0b00,
        0b01,
        0b10
      },{
        0b00,
        0b10,
        0b11
      },{
        0b01,
        0b10,
        0b11
      }
    };
    
    return simplices;
  } else if (dimension == 3) {
    std::vector<std::vector<int>> simplices {
      // Eight "right angle" corners, i.e. a corner and then one step into each direction
      {
        0b000, // anchor
        0b001, // change x1
        0b010, // change x2
        0b100  // change x3
      },{
        0b001, // ...
        0b000,
        0b011,
        0b101
      },{
        0b010,
        0b011,
        0b000,
        0b110
      },{
        0b011,
        0b010,
        0b001,
        0b111
      },{
        0b100,
        0b101,
        0b110,
        0b000
      },{
        0b101,
        0b100,
        0b111,
        0b001
      },{
        0b110,
        0b111,
        0b100,
        0b010
      },{
        0b111,
        0b110,
        0b101,
        0b011
      },{
      // Two simplices that cover the "middle" part of the cube
        0b000,
        0b110,
        0b101,
        0b011
      },{
        0b111,
        0b001,
        0b010,
        0b100
      }
    };
    
    return simplices;
  } else {
    warning("Cannot generate for dimensions other than 2 or 3");
    
    std::vector<std::vector<int>> simplices(0);
    return simplices;
  }
}

std::vector<IntegerVector> getCubeCorners(int dimension) {
  if (dimension == 2) {
    std::vector<IntegerVector> corners = {
      {0, 0},
      {0, 1},
      {1, 0},
      {1, 1}
    };
    
    return corners;
  } else if (dimension == 3) {
    std::vector<IntegerVector> corners = {
      {0, 0, 0},
      {0, 0, 1},
      {0, 1, 0},
      {0, 1, 1},
      {1, 0, 0},
      {1, 0, 1},
      {1, 1, 0},
      {1, 1, 1}
    };
    
    return corners;
  } else {
    warning("Cannot generate for dimensions other than 2 or 3");
    
    std::vector<IntegerVector> corners(0);
    return corners;
  }
}

// [[Rcpp::export]]
List getCriticalPointsCellCPP(NumericMatrix moGradMat, List gradMatList, NumericVector div, IntegerVector locallyNondominated, IntegerVector dims, bool sinks_only) {
  int p = gradMatList.length();
  int d = dims.size();
  int n = div.length();

  std::vector<NumericMatrix> gradMat(p);

  for (int i = 0; i < p; i++) {
    gradMat[i] = as<NumericMatrix>(gradMatList(i));
  }

  std::set<int> locally_nondmominated(locallyNondominated.begin(), locallyNondominated.end());
  
  std::vector<std::vector<int>> simplices = getCellSimplices(d);
  std::vector<IntegerVector> cubeCorners = getCubeCorners(d);
  
  std::unordered_set<int> sinks;
  std::unordered_set<int> sources;
  std::unordered_set<int> saddles;

  std::vector<int> sink_count(n, 0);
  std::vector<int> source_count(n, 0);
  std::vector<int> saddle_count(n, 0);

  // anchor index of the current cell
  IntegerVector anchorIndex;
  IntegerVector neighbourIndex;
  
  std::vector<IntegerVector> cornerIndices(d + 1, d);
  IntegerVector cornerIDs(d + 1);
  
  std::vector<NumericVector> allVectors(p * (d + 1));
  
  std::vector<IntegerVector> cubeIndices(pow(2, d));
  IntegerVector cubeIDs(pow(2, d));
    
  std::vector<NumericVector> cubeVectors(p * pow(2, d));

  // helpers
  NumericVector vec;
  NumericVector compVec;

  for (int id = 1; id <= n; id++) {
    // Loop over cell indices
    // Each anchor index refers to the hypercube defined by itself and +1 into each dimension
    if (id % 1000 == 0) {
      cout << id << '\r';
    }

    anchorIndex = convertCellID2IndicesCPP(id, dims);
    
    if (is_true(any(anchorIndex == dims))) {
      // Hypercube not valid
      continue;
    }
    
    int not_sink_count = 0;
    
    for (int i = 0; i < cubeCorners.size(); i++) {
      cubeIndices[i] = anchorIndex + cubeCorners[i];
      cubeIDs[i] = convertIndices2CellIDCPP(cubeIndices[i], dims);
      
      if (div(cubeIDs[i] - 1) > 0) {
        not_sink_count++;
      }
      
      for (int fID = 0; fID < p; fID++) {
        cubeVectors[p * i + fID] = gradMat[fID](cubeIDs[i] - 1,_);
      }
    }

    if (sinks_only && not_sink_count > (pow(2, d) - (d + 1))) {
      // Skip evaluation if we only want sinks and there cannot be one in this cube
      continue;
    }
    
    int cube_critical = isCritical(cubeVectors);

    if (cube_critical == -1) {
      // If whole cube is not critical, none of its simplices can be
      continue;
    }

    // For each simplex contained in cube
    for (std::vector<int> simplex : simplices) {
      // Collect all indices for nodes in simplex
      bool not_sink = false;
      
      for (int i = 0; i < simplex.size(); i++) {
        cornerIndices[i] = cubeIndices[simplex[i]];
        cornerIDs[i] = cubeIDs[simplex[i]];
        
        if (div(cornerIDs[i] - 1) > 0) {
          not_sink = true;
        }
        
        for (int fID = 0; fID < p; fID++) {
          allVectors[p * i + fID] = gradMat[fID](cornerIDs[i] - 1,_);
        }
      }
      
      if (sinks_only && not_sink) {
        // Skip evaluation if we only want sinks and this cannot be one
        continue;
      }

      bool crit = false;
      
      // Check criticality of MO gradients
      
      int mo_critical = isCritical(allVectors);
      
      if (mo_critical == 1) {
        crit = true;
      } else if (mo_critical == 0) {
        // In the edge case, only critical, if a corner point is locally nondominated
        for (int cornerID: cornerIDs) {
          if (locally_nondmominated.find(cornerID) != locally_nondmominated.end()) {
            crit = true;
          }
        }
      }
      
      // Critical if all point have zero MOG (degenerate critical point)
      
      if (!crit) {
        crit = true;
        
        for (int cornerID : cornerIDs) {
          if (is_true(any(moGradMat(cornerID-1,_) != 0))) {
            crit = false;
            break;
          }
        }
      }

      if (crit) {
        bool pos = false;
        bool neg = false;
        bool zero = false;
        
        for (int i : cornerIDs) {
          if (div(i-1) < 0) {
            neg = true;
          } else if (div(i-1) > 0) {
            pos = true;
          } else {
            zero = true;
          }
        }

        for (int i : cornerIDs) {
          if (pos && !neg) {
            source_count[i-1]++;
          } else if ((neg && !pos) || (zero && !neg && !pos)) {
            // neg && !pos --> regular critical
            // all zero --> degenerate critical
            sink_count[i-1]++;
          } else {
            saddle_count[i-1]++;
          }
        }
      }
    }
  }

  for (int id = 1; id <= n; id++) {
    if (sink_count[id-1] > 0) {
      sinks.insert(id);
    }

    if (source_count[id-1] > 0) {
      sources.insert(id);
    }

    if (saddle_count[id-1] > 0) {
      saddles.insert(id);
    }
  }

  return List::create(Named("sinks")=wrap(sinks), _["sources"]=wrap(sources), _["saddles"]=wrap(saddles));
}

// [[Rcpp::export]]
IntegerVector connectedComponentsGrid(IntegerVector ids, IntegerVector dims) {
  IntegerVector ccs(ids.length(), 0);
  std::unordered_set<int> ids_set(ids.begin(), ids.end());
  std::unordered_map<int,int> id_to_index;

  int d = dims.length();

  for (int i = 0; i < ids.length(); i++) {
    id_to_index[ids(i)] = i;
  }

  int cc = 1;

  IntegerMatrix neighbourhood = getNeighbourhood(d, true);

  for (int start_id : ids) {
    if (ccs(id_to_index[start_id]) == 0) {
      // not yet assigned
      std::deque<int> to_evaluate;
      std::unordered_set<int> seen;

      to_evaluate.push_back(start_id);
      seen.insert(start_id);

      while (!to_evaluate.empty()) {
        int id = to_evaluate.front();
        to_evaluate.pop_front();

        ccs(id_to_index[id]) = cc;

        // add each neighbour to queue, if in ids

        for (int r = 0; r < neighbourhood.nrow(); r++) {
          IntegerVector gridIndex = convertCellID2IndicesCPP(id, dims);
          gridIndex += neighbourhood(r,_);

          if (is_true(any(gridIndex > dims)) || is_true(any(gridIndex < 1))) {
            // invalid grid point
            continue;
          }

          if (is_true(all(neighbourhood(r,_) == 0))) {
            // not interesting
            continue;
          }

          int neighbour_id = convertIndices2CellIDCPP(gridIndex, dims);

          if ((ids_set.find(neighbour_id) != ids_set.end()) && (seen.find(neighbour_id) == seen.end())) {
            to_evaluate.push_back(neighbour_id);
            seen.insert(neighbour_id);
          }
        }
      }

      cc++;
    }

  }

  return ccs;
}

// [[Rcpp::export]]
List integrateVectorField(NumericMatrix gradMat, IntegerVector dims, IntegerVector sinks) {
  std::unordered_set<int> stop(sinks.begin(), sinks.end());

  NumericVector currentPoint;
  NumericVector currentPointCell;
  IntegerVector currentIndex;
  int currentID;
  int nextID;
  int previousID = -1;

  int d = dims.size();
  int n = gradMat.nrow();

  // helpers
  double scaling;
  double requiredScaling;
  double target;

  std::vector<double> height(n);
  std::vector<int> last_visited(n);
  std::vector<double> gradientLengths(n);
  std::vector< std::vector<double> > gradients(n);

  for (int id = 1; id <= n; id++) {
    gradientLengths[id-1] = computeVectorLengthCPP(gradMat(id-1,_));

    if (stop.find(id) != stop.end()) {
      height[id-1] = gradientLengths[id-1] / 2;
      last_visited[id-1] = id;
    }
  }

  for (int id = 1; id <= n; id++) {
    if (id % 1000 == 0) {
      cout << '\r' << id;
    }

    std::unordered_set<int> seen;

    if (stop.find(id) != stop.end()) continue;

    currentID = id;

    currentIndex = convertCellID2IndicesCPP(id, dims);
    currentPoint = as<NumericVector>(currentIndex);

    while (true) {
      scaling = std::numeric_limits<double>::max();

      for (int d_iter = 0; d_iter < d; d_iter++) {
        double grad_d_iter = gradMat(currentID-1,d_iter);

        if (grad_d_iter > 0) {
          // check collision with cell + 1 in this dimension
          target = currentIndex(d_iter) + 0.51;
        } else if (grad_d_iter < 0) {
          // check collision with cell - 1 in this dimension
          target = currentIndex(d_iter) - 0.51;
        } else {
          // grad_d_iter == 0
          continue;
        }

        requiredScaling = (target - currentPoint(d_iter)) / grad_d_iter;

        scaling = min(requiredScaling, scaling);
      }

      if (scaling == std::numeric_limits<double>::max()) {
        break;
      }

      currentPoint = currentPoint + scaling * gradMat(currentID-1,_);
      height[id - 1] += scaling * gradientLengths[currentID - 1] * gradientLengths[currentID - 1];

      currentPointCell = round(currentPoint, 0);
      currentIndex = as<IntegerVector>(currentPointCell);

      if (is_true(any(currentIndex > dims)) || is_true(any(currentIndex < 1))) {
        break;
      }

      if (d == 2) {
        nextID = currentIndex(0) + (currentIndex(1) - 1) * dims[0];
      } else {
        nextID = convertIndices2CellIDCPP(currentIndex, dims);
      }

      if (nextID != previousID && seen.find(nextID) != seen.end()) {
        break;
      }

      seen.insert(currentID);
      previousID = currentID;
      currentID = nextID;

      last_visited[id-1] = currentID;

      if (stop.find(currentID) != stop.end()) {
        break;
      }
    }
  }

  return(List::create(_["height"]=height, _["last.visited"]=last_visited));
}

// [[Rcpp::export]]
IntegerVector locallyNondominatedCPP(NumericMatrix fnMat, IntegerVector dims, bool includeDiagonals) {
  int n = fnMat.nrow();
  int d = dims.size();

  std::vector<int> locallyNondominated;

  IntegerMatrix neighbourhood = getNeighbourhood(d, includeDiagonals);

  int neighbours = neighbourhood.nrow();

  // helper
  IntegerVector indices;
  IntegerVector neighbourIndices;
  bool dominated;
  int neighbourID;

  for (int id = 1; id <= n; id++) {
    indices = convertCellID2IndicesCPP(id, dims);
    dominated = false;

    for (int r = 0; r < neighbours; r++) {
      neighbourIndices = indices + neighbourhood(r, _);
      if (is_true(any(neighbourIndices < 1)) ||
          is_true(any(neighbourIndices > dims)) ||
          is_true(all(neighbourIndices == indices))) {
        continue;
      }

      neighbourID = convertIndices2CellIDCPP(neighbourIndices, dims);

      if (is_true(any(fnMat(neighbourID - 1, _) < fnMat(id - 1, _))) &&
          is_true(all(fnMat(neighbourID - 1, _) <= fnMat(id - 1, _)))) {
        dominated = true;
        break;
      }
    }

    if (!dominated) {
      locallyNondominated.push_back(id);
    }
  }

  return wrap(locallyNondominated);
}

// [[Rcpp::export]]
IntegerVector changeOfSignCPP(NumericVector fnVec, IntegerVector dims, bool includeDiagonals) {
  int n = fnVec.size();
  int d = dims.size();

  IntegerVector changeOfSignIndices;

  IntegerMatrix deltas = getNeighbourhood(d, includeDiagonals);

  int nDeltas = deltas.nrow();

  // helper
  IntegerVector indices;
  IntegerVector neighbourIndices;
  bool signChange;
  int neighbourID;

  for (int id = 1; id <= n; id++) {
    if (fnVec(id - 1) == 0) {
      changeOfSignIndices.push_back(id);
      continue;
    }

    indices = convertCellID2IndicesCPP(id, dims);
    signChange = false;

    for (int r = 0; r < nDeltas; r++) {
      neighbourIndices = indices + deltas(r, _);
      if (is_true(any(neighbourIndices < 1)) ||
          is_true(any(neighbourIndices > dims))) {
        // invalid neighbours
        continue;
      }

      neighbourID = convertIndices2CellIDCPP(neighbourIndices, dims);

      if (fnVec(id - 1) * fnVec(neighbourID - 1) < 0) {
        // some neighbour has a different sign!
        signChange = true;
        break;
      }
    }

    if (signChange) {
      changeOfSignIndices.push_back(id);
    }
  }

  return changeOfSignIndices;
}

// [[Rcpp::export]]
IntegerVector changeOfBasin(IntegerVector basins, IntegerVector dims) {
  int n = basins.size();
  int d = dims.size();

  std::vector<int> ridges;

  IntegerMatrix neighbourhood = getNeighbourhood(d, false);

  int neighbours_amt = neighbourhood.nrow();

  // helper
  IntegerVector indices;
  IntegerVector neighbourIndices;
  bool basin_change;
  int neighbourID;

  for (int id = 1; id <= n; id++) {
    indices = convertCellID2IndicesCPP(id, dims);
    basin_change = false;

    for (int r = 0; r < neighbours_amt; r++) {
      neighbourIndices = indices + neighbourhood(r, _);
      if (is_true(any(neighbourIndices < 1)) ||
          is_true(any(neighbourIndices > dims))) {
        // invalid neighbour
        continue;
      }

      neighbourID = convertIndices2CellIDCPP(neighbourIndices, dims);

      if (basins(id - 1) != basins(neighbourID - 1)) {
        // some neighbour is in a different basin!
        basin_change = true;
        break;
      }
    }

    if (basin_change) {
      ridges.push_back(id);
    }
  }

  return wrap(ridges);
}

// [[Rcpp::export]]
NumericMatrix gridBasedGradientCPP(NumericVector fnVec, IntegerVector dims, NumericVector stepSizes, double precNorm, double precAngle) {
  int n = fnVec.size();
  int d = dims.size();
  // FIXME does not work properly!
  bool scaling = (stepSizes != nullptr);
  NumericMatrix gradMat(n, d);

  // helper variables
  double diff;
  IntegerVector indices;
  IntegerVector iPlus;
  IntegerVector iMinus;
  double f;
  double fPlus;
  double fMinus;

  for (int id = 1; id <= n; id++) {
    indices = convertCellID2IndicesCPP(id, dims);

    for (int dim = 0; dim < d; dim++) {
      // vector access is zero-indexed!
      if (indices(dim) == 1) {
        iPlus = clone(indices);
        iPlus(dim)++;

        fPlus = fnVec[convertIndices2CellIDCPP(iPlus, dims) - 1];
        f = fnVec[convertIndices2CellIDCPP(indices, dims) - 1];

        diff = fPlus - f;
      } else if (indices(dim) == dims(dim)) {
        iMinus = clone(indices);
        iMinus(dim)--;

        fMinus = fnVec[convertIndices2CellIDCPP(iMinus, dims) - 1];
        f = fnVec[convertIndices2CellIDCPP(indices, dims) - 1];

        diff = f - fMinus;
      } else {
        iPlus = clone(indices);
        iPlus(dim)++;
        iMinus = clone(indices);
        iMinus(dim)--;

        fPlus = fnVec[convertIndices2CellIDCPP(iPlus, dims) - 1];
        fMinus = fnVec[convertIndices2CellIDCPP(iMinus, dims) - 1];

        diff = (fPlus - fMinus) / 2;
      }

      if (scaling) {
        diff = diff / stepSizes(dim);
      }

      gradMat(id-1, dim) = diff;
    }
  }

  return gradMat;
}

// [[Rcpp::export]]
List cumulateGradientsCPP(NumericMatrix centers, NumericMatrix gradients, IntegerVector stopCells, double precVectorLength, double precNorm, bool fixDiagonals, bool cumulateGradientLength) {
  int d = centers.ncol();                   // dimensionality of search space
  int n = centers.nrow();                   // number of grid points

  // number of grid points per (search space) dimension
  IntegerVector dims(d, 0);
  for (int j = 0; j < d; j++) {
    NumericVector ctr = centers(_,j);
    ctr = unique(ctr);
    dims[j] = ctr.size();
  }

  NumericVector gradientLengths(n);         // length vector for the multi-objective gradients
  IntegerVector cellPointer(n, -999);       // vector, which indicates per cell the successor cell
  NumericVector gradFieldVector(n, -999.0); // the "final" result vector containing the cell height
  IntegerVector last_visited(n, -1);        // vector containing the cell that integration was stopped at
  LogicalVector visited = rep(false, n);    // which cells have already been "visited" / "processed"

  // helper variables
  double vectorLength = -1.0;
  NumericVector currentGradients;
  int visitCounter = 0;                     // counter, enabling early stop of algorithm once all cells are processed
  IntegerVector currentCell(d);             // row and column index of the current cell
  IntegerVector nextCell(d);                // row and column index of the next (= successor) cell
  int nextCellID;

  for (int i = 0; i < n; i++) {
    // iterate over all cells, extract their gradients and store the gradient lengths
    currentGradients = gradients(i,_);
    vectorLength = computeVectorLengthCPP(currentGradients);
    gradientLengths[i] = vectorLength;
    if (vectorLength < precVectorLength) {
      // if the multi-objective gradient is "short enough", the current cell is considered
      // to be locally efficient (--> move on to the next cell)
      visited[i] = true;
      visitCounter++;
      gradFieldVector[i] = vectorLength;
    } else {
      // nextCell = findNextCellHighDimensionCPP(currentGradients);
      // currentCell = convertHighDimensionalCellID2IndicesCPP(i + 1, dims);
      nextCell = findNextCellCPP(currentGradients);
      currentCell = convertCellID2IndicesCPP(i + 1, dims);
      nextCell += currentCell;

      if (is_true(any(nextCell > dims)) | is_true(any(nextCell < 1))) {
        // if the next cell is located outside the boundaries, one can not move further
        // into the direction of its gradient
        visited[i] = true;
        visitCounter++;
        gradFieldVector[i] = vectorLength;
        last_visited[i] = i+1; // id = i+1 here
      } else {
        // otherwise, store the successor cell
        // nextCellID = convertHighDimensionalIndices2CellIDCPP(nextCell, dims);
        nextCellID = convertIndices2CellIDCPP(nextCell, dims);
        cellPointer[i] = nextCellID - 1;
      }
    }
  }

  // set 0 height and visited for each stop cell
  for (int id : stopCells) {
    visitCounter++;
    visited[id-1] = true;
    last_visited[id-1] = id;
    gradFieldVector(id-1) = gradientLengths(id-1) / 2;
  }

  /* below, the actual cumulative part begins */

  // helper variables
  int currentCellID;
  int pathLength;
  IntegerVector path(n, -999);
  for (int i = 0; i < n; i++) {
    // iterate over all cells
    if (visitCounter == n) {
      // leave the for-loop if all cells have already been visited
      break;
    }
    if (!visited[i]) {
      // if the i-th cell has not yet been visited, follow the path
      // from that cell towards a previously "visited" cell (either
      // from an earlier path, a boundary point pointing out of bounds,
      // or a local efficient point)

      // initialize the path with cell i
      currentCellID = i;
      path[0] = currentCellID;
      pathLength = 1;
      visited[currentCellID] = true;
      nextCellID = cellPointer[currentCellID];
      while (!visited[nextCellID]) {
        // as long as the succeeding cell has not been visited, append
        // such a cell to the current path
        currentCellID = nextCellID;
        path[pathLength] = currentCellID;
        visited[currentCellID] = true;
        visitCounter++;
        nextCellID = cellPointer[currentCellID];
        pathLength++;
      }

      double a = 0.0;
      if (gradFieldVector[nextCellID] > -999) {
        // if the path stopped, because its next cell was already
        // part of an earlier loop, add the value of the next
        // cell to the cumulated sum (a) -- and ensure that
        // this cell is not counted twice
        a = gradFieldVector[nextCellID];
        visitCounter--;
      }
      if (pathLength > 0) {
        // iterate backwards over the path and add the cumulated
        // sum of gradient lengths to the cells of the path
        for (int j = pathLength; j > 0; j--) {
          int index = path[j - 1];
          double length;

          if (cumulateGradientLength) {
            length = gradientLengths[index];
          } else {
            length = 1.0;
          }

          if(fixDiagonals) {
            NumericVector delta = as<NumericVector>(findNextCellCPP(gradients(index,_)));
            double stepSize = computeVectorLengthCPP(delta);
            a += length * stepSize;
          } else {
            a += length;
          }

          gradFieldVector[index] = a;
          if (last_visited[nextCellID] != -1) {
            last_visited[index] = last_visited[nextCellID];
          } else {
            last_visited[index] = currentCellID + 1;
          }

        }
      }
    }
  }

  return List::create(_["height"]=gradFieldVector, _["last.visited"]=last_visited);
}

// [[Rcpp::export]]
NumericVector getBiObjGradientCPP(NumericVector g1, NumericVector g2, double precNorm, double precAngle) {
  // Multiobjective gradient for two objectives
  int len = g1.length();
  NumericVector zeros (len);
  zeros.attr("dim") = R_NilValue;

  g1 = normalizeVectorCPP(g1, precNorm);
  if (is_true(all(g1 == 0))) {
    // if the gradient of fn1 is zero, this has to be a local efficient point
    return(zeros);
  }

  g2 = normalizeVectorCPP(g2, precNorm);
  if (is_true(all(g2 == 0))) {
    // if the gradient of fn2 is zero, this has to be a local efficient point
    return(zeros);
  }

  double angle1 = computeAngleCPP(g1, g2, precNorm);

  if (abs(180 - angle1) < precAngle) {
    // if the angle between both gradients is (approximately) 180 degree,
    // this has to be a local efficient point
    return(zeros);
  }

  return(0.5 * (g1 + g2));
}

// [[Rcpp::export]]
NumericVector getTriObjGradientCPP(NumericVector g1, NumericVector g2, NumericVector g3, double precNorm, double precAngle) {
  // Multiobjective gradient for three objectives
  int len = g1.length();
  NumericVector zeros (len);

  g1 = normalizeVectorCPP(g1, precNorm);
  g2 = normalizeVectorCPP(g2, precNorm);
  g3 = normalizeVectorCPP(g3, precNorm);

  if (is_true(all(g1 == 0)) ||
      is_true(all(g2 == 0)) ||
      is_true(all(g3 == 0))) {
    // if the gradient of any objective is zero, this has to be a local efficient point
    return(zeros);
  }

  double angle1 = computeAngleCPP(g1, g2, precNorm);
  double angle2 = computeAngleCPP(g1, g3, precNorm);
  double angle3 = computeAngleCPP(g2, g3, precNorm);

  if (abs(180 - angle1) < precAngle ||
      abs(180 - angle2) < precAngle ||
      abs(180 - angle3) < precAngle) {
    // if the angle between any two gradients is (approximately) 180 degree,
    // this has to be a local efficient point
    return(zeros);
  }

  if (abs(angle1 + angle2 + angle3 - 360) < precAngle) {
    // if all gradients show in "opposite" directions, this has to be a local effient point
    return(zeros);
  }

  std::vector<NumericVector> vectors;
  vectors.push_back(g1);
  vectors.push_back(g2);
  vectors.push_back(g3);

  if (len == 2 && isCritical(vectors) == 1) {
    return(zeros);
  }

  if (len >= 3) {
    NumericVector n = normalizeVectorCPP(compute3DcrossProductCPP(g1-g2,g1-g3), 0);
    NumericVector moNorm = n * sum(n * g1); // n * (n dot g1) = n * (distance to plane created by g1,g2,g3)

    double normAngle1 = computeAngleCPP(g1 - moNorm, g2 - moNorm, precNorm);
    double normAngle2 = computeAngleCPP(g1 - moNorm, g3 - moNorm, precNorm);
    double normAngle3 = computeAngleCPP(g2 - moNorm, g3 - moNorm, precNorm);

    if (abs(normAngle1 + normAngle2 + normAngle3 - 360) <= precAngle) {
      // norm is a convex combination of
      // single-objective gradients
      return(moNorm);
    } else {
      // if the norm vector is not a convex combination,
      // the MOG is pointing to one of the boundaries
      // --> follow widest angle
      double maxAngle = max(NumericVector::create(angle1, angle2, angle3));
      if (angle1 == maxAngle) {
        return(0.5 * (g1 + g2));
      } else if (angle2 == maxAngle) {
        return(0.5 * (g1 + g3));
      } else {
        return(0.5 * (g2 + g3));
      }
    }
  } else {
    double maxAngle = max(NumericVector::create(angle1, angle2, angle3));
    if (angle1 == maxAngle) {
      return(0.5 * (g1 + g2));
    } else if (angle2 == maxAngle) {
      return(0.5 * (g1 + g3));
    } else {
      return(0.5 * (g2 + g3));
    }
  }
}

// [[Rcpp::export]]
NumericMatrix getBiObjGradientGridCPP(NumericMatrix gradMat1, NumericMatrix gradMat2, double precNorm, double precAngle) {
  int n = gradMat1.rows();
  int d = gradMat1.cols();
  NumericMatrix moGradMat(n, d);
  for (int i = 0; i < n; i++) {
    moGradMat(i,_) = getBiObjGradientCPP(gradMat1(i,_), gradMat2(i,_), precNorm, precAngle);
  }
  return moGradMat;
}

// [[Rcpp::export]]
NumericMatrix getTriObjGradientGridCPP(NumericMatrix gradMat1, NumericMatrix gradMat2, NumericMatrix gradMat3, double precNorm, double precAngle) {
  int n = gradMat1.rows();
  int d = gradMat1.cols();
  NumericMatrix moGradMat(n, d);
  for (int i = 0; i < n; i++) {
    moGradMat(i,_) = getTriObjGradientCPP(gradMat1(i,_), gradMat2(i,_), gradMat3(i,_), precNorm, precAngle);
  }
  return moGradMat;
}

// [[Rcpp::export]]
NumericVector calculateMaxDisplayHeightCPP(NumericVector heights, IntegerVector dims, bool includeDiagonals) {
  // points are "dominated" by their neighbours from some point on
  // calculate this point here
  int n = heights.length();
  int d = dims.length();
  NumericVector maxHeights(n);
  IntegerMatrix deltas = getNeighbourhood(d, includeDiagonals);
  int nDeltas = deltas.nrow();
  IntegerVector indices;
  IntegerVector neighbourIndices;
  int neighbourID;
  double maxNeighbour;
  for (int id = 1; id <= n; id++) {
    indices = convertCellID2IndicesCPP(id, dims);
    maxNeighbour = 0.0;
    for (int r = 0; r < nDeltas; r++) {
      neighbourIndices = indices + deltas(r, _);
      if (is_true(any(neighbourIndices < 1)) ||
          is_true(any(neighbourIndices > dims))) {
        // keep Infinity as maxHeight when at boundaries
        maxNeighbour = std::numeric_limits<double>::infinity();
        break;
      }
      neighbourID = convertIndices2CellIDCPP(neighbourIndices, dims);
      maxNeighbour = std::max(maxNeighbour, heights(neighbourID-1));
    }
    maxHeights(id-1) = maxNeighbour;
  }
  return maxHeights;
}

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>
#define INF 1e5 + 1
int32_t Sign(int64_t val) {
  int32_t sign = 0;
  if (val < 0) {
    sign = -1;
  }
  if (val > 0) {
    sign = 1;
  }
  return sign;
}
namespace Geometry {
class Vector;

class IShape;
class Point;
class Segment;
class Ray;
class Line;
class Circle;
class Polygon;
class IShape {
 public:
  virtual IShape& Move(const Vector& vec) = 0;
  [[nodiscard]] virtual bool ContainsPoint(const Point& point) const = 0;
  [[nodiscard]] virtual bool CrossesSegment(const Segment& seg) const = 0;
  [[nodiscard]] virtual IShape* Clone() const = 0;
  virtual std::string ToString() = 0;
  virtual ~IShape() = default;
};

class Point : public IShape {
 public:
  IShape& Move(const Vector& vec) override;
  [[nodiscard]] bool ContainsPoint(const Point& point) const override;
  [[nodiscard]] bool CrossesSegment(const Segment& seg) const override;
  [[nodiscard]] IShape* Clone() const override;
  std::string ToString() override;
  explicit Point(int64_t x = 0, int64_t y = 0);
  ~Point() override = default;
  int64_t x;
  int64_t y;
};

class Vector {
 public:
  Vector();
  Vector(const Point& begin, const Point& end);
  Vector operator+(const Vector& second) const;
  int64_t operator*(const Vector& second) const;
  int64_t operator^(const Vector& second) const;
  [[nodiscard]] double Length() const;
  Vector operator-() const;
  [[nodiscard]] Point ReturnA() const;
  ~Vector() = default;

 private:
  Point a_;
};

Vector operator-(const Point& first, const Point& second) {
  Point tmp(first.x, first.y);
  tmp.x -= second.x;
  tmp.y -= second.y;
  Vector ret = Vector(Point(0, 0), tmp);
  return ret;
}

class Line : public IShape {
 public:
  Line(const Point& a, const Point& b);
  [[nodiscard]] bool CheckParallel(const Line& line) const;
  [[nodiscard]] Vector ReturnDirectionVector() const;
  IShape& Move(const Vector& vec) override;
  [[nodiscard]] bool ContainsPoint(const Point& point) const override;
  [[nodiscard]] bool CrossesSegment(const Segment& seg) const override;
  [[nodiscard]] IShape* Clone() const override;
  std::string ToString() override;
  [[nodiscard]] std::pair<double, double> FindIntersectionPoint(
      const Line& line) const;
  [[nodiscard]] double FindDist(const Line& line) const;
  [[nodiscard]] double FindDistWithPoint(const Point& point) const;
  ~Line() override = default;

 private:
  int64_t a_ = 0;
  int64_t b_ = 0;
  int64_t c_ = 0;
  Point first_ = Point(0, 0);
  Point second_ = Point(0, 0);
};

class Ray : public IShape {
 public:
  Ray(const Point& a, const Point& b);
  double FindDistWithPoint(Point& c) const;
  IShape& Move(const Vector& vec) override;
  [[nodiscard]] bool ContainsPoint(const Point& c) const override;
  [[nodiscard]] bool CrossesSegment(const Segment& seg) const override;
  [[nodiscard]] IShape* Clone() const override;
  std::string ToString() override;
  ~Ray() override = default;

 private:
  Point a_ = Point(0, 0);
  Point b_ = Point(0, 0);
};

class Segment : public IShape {
 public:
  Segment(const Point& a, const Point& b);
  [[nodiscard]] Point ReturnBegin() const;
  [[nodiscard]] Point ReturnEnd() const;
  [[nodiscard]] double FindDistWithPoint(const Point& c) const;
  IShape& Move(const Vector& vec) override;
  [[nodiscard]] bool ContainsPoint(const Point& c) const override;
  [[nodiscard]] bool CrossesSegment(const Segment& seg) const override;
  [[nodiscard]] IShape* Clone() const override;
  std::string ToString() override;
  ~Segment() override = default;

 private:
  Point a_ = Point(0, 0);
  Point b_ = Point(0, 0);
};

class Polygon : public IShape {
 public:
  explicit Polygon(const std::vector<Point>& vertexes);
  bool IsConvex();
  IShape& Move(const Vector& vec) override;
  [[nodiscard]] bool ContainsPoint(const Point& c) const override;
  [[nodiscard]] bool CrossesSegment(const Segment& seg) const override;
  [[nodiscard]] IShape* Clone() const override;
  std::string ToString() override;
  ~Polygon() override = default;

 private:
  std::vector<Point> vertexes_;
};

class Circle : public IShape {
 public:
  explicit Circle(const Point& centre, int64_t radius);
  IShape& Move(const Vector& vec) override;
  [[nodiscard]] bool ContainsPoint(const Point& c) const override;
  [[nodiscard]] bool CrossesSegment(const Segment& seg) const override;
  [[nodiscard]] IShape* Clone() const override;
  std::string ToString() override;
  ~Circle() override = default;

 private:
  Point centre_ = Point(0, 0);
  int64_t radius_ = 0;
};

Point::Point(int64_t x, int64_t y) {
  this->x = x;
  this->y = y;
}

IShape& Point::Move(const Vector& vec) {
  x += vec.ReturnA().x;
  y += vec.ReturnA().y;
  return *this;
}

bool Point::ContainsPoint(const Point& point) const {
  bool contains = false;
  if (x == point.x && y == point.y) {
    contains = true;
  }
  return contains;
}

IShape* Point::Clone() const {
  auto* cp = new Point(x, y);
  return cp;
}

bool Point::CrossesSegment(const Segment& seg) const {
  Point tmp(x, y);
  bool crosses = false;
  if (seg.ContainsPoint(tmp)) {
    crosses = true;
  }
  return crosses;
}

std::string Point::ToString() {
  std::string point_str =
      "Point(" + std::to_string(x) + ", " + std::to_string(y) + ")";
  return point_str;
}

Vector::Vector() = default;

Vector::Vector(const Point& begin, const Point& end) {
  a_.x = end.x - begin.x;
  a_.y = end.y - begin.y;
}
int64_t Vector::operator*(const Vector& second) const {
  int64_t scalar_mul = a_.x * second.a_.x + a_.y * second.a_.y;
  return scalar_mul;
}
int64_t Vector::operator^(const Vector& second) const {
  int64_t vec_mul = a_.x * second.a_.y - a_.y * second.a_.x;
  return vec_mul;
}
Vector Vector::operator-() const {
  Point begin = Point(0, 0);
  Point new_end = Point(-a_.x, -a_.y);
  Vector ret = Vector(begin, new_end);
  return ret;
}
Vector Vector::operator+(const Vector& second) const {
  Point start = Point(0, 0);
  Point end = Point(a_.x + second.a_.x, a_.y + second.a_.y);
  Vector sum = Vector(start, end);
  return sum;
}
double Vector::Length() const {
  return sqrt(static_cast<double>(a_.x * a_.x + a_.y * a_.y));
}

Point Vector::ReturnA() const { return a_; }

Vector Line::ReturnDirectionVector() const {
  Point start = Point(0, 0);
  Point end = Point(b_, -a_);
  return {start, end};
}
std::pair<double, double> Line::FindIntersectionPoint(const Line& line) const {
  int64_t det = a_ * line.b_ - b_ * line.a_;
  int64_t det_x = b_ * line.c_ - c_ * line.b_;
  int64_t det_y = -a_ * line.c_ + c_ * line.a_;
  double x = static_cast<double>(det_x) / static_cast<double>(det);
  double y = static_cast<double>(det_y) / static_cast<double>(det);
  return std::make_pair(x, y);
}
double Line::FindDist(const Line& line) const {
  double proportion = 0;
  if (a_ != 0) {
    proportion = static_cast<double>(line.a_) / static_cast<double>(a_);
  } else {
    proportion = static_cast<double>(line.b_) / static_cast<double>(b_);
  }
  double c1 = static_cast<double>(c_) * proportion;
  return std::abs(c1 - static_cast<double>(line.c_)) /
         sqrt(static_cast<double>(a_ * a_ + b_ * b_));
}

Line::Line(const Point& a, const Point& b) {
  this->a_ = b.y - a.y;
  this->b_ = a.x - b.x;
  this->c_ = a.y * (b.x - a.x) - a.x * (b.y - a.y);
  first_ = a;
  second_ = b;
}

bool Line::ContainsPoint(const Point& point) const {
  bool contains = false;
  int64_t res = a_ * point.x + b_ * point.y + c_;
  if (res == 0) {
    contains = true;
  }
  return contains;
}

double Line::FindDistWithPoint(const Point& point) const {
  return static_cast<double>(
      std::abs(static_cast<double>(a_) * static_cast<double>(point.x) +
               static_cast<double>(b_) * static_cast<double>(point.y) +
               static_cast<double>(c_)) /
      sqrt(static_cast<double>(a_ * a_) + static_cast<double>(b_ * b_)));
}

bool Line::CheckParallel(const Line& line) const {
  bool parallel = false;
  Vector dir_first = ReturnDirectionVector();
  Vector dir_second = line.ReturnDirectionVector();
  int64_t vec_mul = dir_second ^ dir_first;
  if (vec_mul == 0) {
    parallel = true;
  }
  return parallel;
}
IShape& Line::Move(const Vector& vec) {
  c_ = c_ - a_ * vec.ReturnA().x - b_ * vec.ReturnA().y;
  first_.x += vec.ReturnA().x;
  first_.y += vec.ReturnA().y;
  second_.x += vec.ReturnA().x;
  second_.y += vec.ReturnA().y;
  return *this;
}
std::string Line::ToString() {
  std::string line_str = "Line(" + std::to_string(a_) + ", " +
                         std::to_string(b_) + ", " + std::to_string(c_) + ")";
  return line_str;
}
IShape* Line::Clone() const {
  auto* cp = new Line(first_, second_);
  return cp;
}
bool Line::CrossesSegment(const Segment& seg) const {
  bool crosses = false;
  Vector ab = Vector(first_, second_);
  Vector ac = Vector(first_, seg.ReturnBegin());
  Vector ad = Vector(first_, seg.ReturnEnd());
  int64_t sign_ab_ac = Sign(ab ^ ac);
  int64_t sign_ab_ad = Sign(ab ^ ad);
  if (sign_ab_ac * sign_ab_ad <= 0) {
    crosses = true;
  }
  return crosses;
}

Ray::Ray(const Point& a, const Point& b) {
  this->a_ = a;
  this->b_ = b;
}

bool Ray::ContainsPoint(const Point& c) const {
  bool contains = false;
  Vector first = Vector(a_, b_);
  Vector second = Vector(a_, c);
  int64_t scal_mul = Sign(first * second);
  if (Sign(first ^ second) == 0 && scal_mul >= 0) {
    contains = true;
  }
  return contains;
}

double Ray::FindDistWithPoint(Point& c) const {
  Vector ab = Vector(a_, b_);
  Vector ac = Vector(a_, c);
  Line tmp = Line(a_, b_);
  int64_t sc_mul = ab * ac;
  double ret = 0;
  if (sc_mul >= 0) {
    ret = tmp.FindDistWithPoint(c);
  } else {
    ret = ac.Length();
  }
  return ret;
}
IShape& Ray::Move(const Vector& vec) {
  a_.x += vec.ReturnA().x;
  a_.y += vec.ReturnA().y;
  b_.x += vec.ReturnA().x;
  b_.y += vec.ReturnA().y;
  return *this;
}
bool Ray::CrossesSegment(const Segment& seg) const {
  bool crosses = false;
  Vector ray_vec = Vector(a_, b_);
  Vector seg_vec = Vector(seg.ReturnBegin(), seg.ReturnEnd());
  Vector from_seg_to_ray = Vector(seg.ReturnBegin(), a_);
  if (Sign(ray_vec ^ seg_vec) == 0) {
    if (ContainsPoint(seg.ReturnBegin()) || ContainsPoint(seg.ReturnEnd())) {
      crosses = true;
    }
    return crosses;
  }
  int64_t rs = ray_vec ^ seg_vec;
  int64_t rsr = ray_vec ^ from_seg_to_ray;
  int64_t srs = seg_vec ^ from_seg_to_ray;
  if (rs < 0) {
    rs *= -1;
    rsr *= -1;
    srs *= -1;
  }
  if (rsr >= 0 && rsr <= rs && srs >= 0) {
    crosses = true;
  }
  return crosses;
}
IShape* Ray::Clone() const {
  auto* cp = new Ray(a_, b_);
  return cp;
}
std::string Ray::ToString() {
  Vector tmp = Vector(a_, b_);
  std::string ray_str = "Ray(" + a_.ToString() + ", " + "Vector(" +
                        std::to_string(tmp.ReturnA().x) + ", " +
                        std::to_string(tmp.ReturnA().y) + "))";
  return ray_str;
}

Segment::Segment(const Point& a, const Point& b) {
  this->a_ = a;
  this->b_ = b;
}

double Segment::FindDistWithPoint(const Point& c) const {
  Vector ac = Vector(a_, c);
  Vector ab = Vector(a_, b_);
  Vector bc = Vector(b_, c);
  Line tmp = Line(a_, b_);
  int64_t sc_first = ac * ab;
  int64_t sc_sec = bc * (-ab);
  double ret = 0;
  if (sc_first >= 0 && sc_sec >= 0) {
    ret = tmp.FindDistWithPoint(c);
  } else {
    ret = std::min(ac.Length(), bc.Length());
  }
  return ret;
}
bool Segment::CrossesSegment(const Segment& seg) const {
  bool crosses = false;
  Vector ab = Vector(a_, b_);
  Vector ac = Vector(a_, seg.a_);
  Vector ad = Vector(a_, seg.b_);
  Vector cd = Vector(seg.a_, seg.b_);
  Vector cb = Vector(seg.a_, b_);
  Vector db = Vector(seg.b_, b_);
  int64_t vec_mul_ab_ac = ab ^ ac;
  int64_t vec_mul_ab_ad = ab ^ ad;
  int64_t vec_mul_cd_ca = -ac ^ cd;
  int64_t vec_mul_cd_cb = cb ^ cd;
  if ((Sign(vec_mul_ab_ac) * Sign(vec_mul_ab_ad) < 0 &&
       Sign(vec_mul_cd_ca) * Sign(vec_mul_cd_cb) < 0) ||
      (Sign(vec_mul_ab_ac) == 0 && Sign(-ac * cb) <= 0) ||
      (Sign(vec_mul_ab_ad) == 0 && Sign(-ad * db) <= 0) ||
      (Sign(vec_mul_cd_ca) == 0 && Sign(ac * ad) <= 0) ||
      (Sign(vec_mul_cd_cb) == 0 && Sign(-cb * -db) <= 0)) {
    crosses = true;
  }
  return crosses;
}
bool Segment::ContainsPoint(const Point& c) const {
  bool contains = false;
  Vector ac = Vector(a_, c);
  Vector ab = Vector(a_, b_);
  int64_t vec_mul = ab ^ ac;
  if (vec_mul == 0) {
    if (c.x <= std::max(a_.x, b_.x) && c.x >= std::min(a_.x, b_.x) &&
        c.y <= std::max(a_.y, b_.y) && c.y >= std::min(a_.y, b_.y)) {
      contains = true;
    }
    return contains;
  }
  return contains;
}
IShape& Segment::Move(const Vector& vec) {
  a_.x += vec.ReturnA().x;
  a_.y += vec.ReturnA().y;
  b_.x += vec.ReturnA().x;
  b_.y += vec.ReturnA().y;
  return *this;
}
IShape* Segment::Clone() const {
  auto* cp = new Segment(a_, b_);
  return cp;
}
std::string Segment::ToString() {
  std::string seg_str = "Segment(" + a_.ToString() + ", " + b_.ToString() + ")";
  return seg_str;
}
Point Segment::ReturnBegin() const { return a_; }
Point Segment::ReturnEnd() const { return b_; }

bool Polygon::IsConvex() {
  bool convex = false;
  int32_t sign = 0;
  uint32_t size = vertexes_.size();
  Vector first = Vector(vertexes_[size - 1], vertexes_[0]);
  Vector second = Vector(vertexes_[0], vertexes_[1]);
  sign = Sign(first ^ second);
  for (uint32_t i = 1; i < size - 1; i++) {
    first = Vector(vertexes_[i - 1], vertexes_[i]);
    second = Vector(vertexes_[i], vertexes_[i + 1]);
    if (Sign(first ^ second) * sign < 0) {
      return convex;
    }
    if (Sign(first ^ second) != 0 && sign == 0) {
      sign = Sign(first ^ second);
    }
  }
  first = Vector(vertexes_[size - 2], vertexes_[size - 1]);
  second = Vector(vertexes_[size - 1], vertexes_[0]);
  if (Sign(first ^ second) * sign < 0) {
    return convex;
  }
  convex = true;
  return convex;
}
Polygon::Polygon(const std::vector<Point>& vertexes) {
  this->vertexes_ = vertexes;
}
IShape& Polygon::Move(const Vector& vec) {
  for (auto& vertex : vertexes_) {
    vertex.x += vec.ReturnA().x;
    vertex.y += vec.ReturnA().y;
  }
  return *this;
}
IShape* Polygon::Clone() const {
  auto* cp = new Polygon(vertexes_);
  return cp;
}
std::string Polygon::ToString() {
  std::string poly_str = "Polygon(";
  for (uint32_t i = 0; i < vertexes_.size() - 1; ++i) {
    poly_str += vertexes_[i].ToString() + ", ";
  }
  poly_str += vertexes_[vertexes_.size() - 1].ToString() + ")";
  return poly_str;
}
bool Polygon::ContainsPoint(const Point& c) const {
  bool contains = true;
  Point begin = Point(c.x, c.y);
  Point end = Point(INF, c.y + 1);
  Segment curr = Segment(vertexes_[vertexes_.size() - 1], vertexes_[0]);
  if (curr.ContainsPoint(c)) {
    return contains;
  }
  for (uint32_t i = 0; i < vertexes_.size() - 1; ++i) {
    curr = Segment(vertexes_[i], vertexes_[i + 1]);
    if (curr.ContainsPoint(c)) {
      return contains;
    }
  }
  Segment tmp = Segment(begin, end);
  int32_t intersects = 0;
  curr = Segment(vertexes_[vertexes_.size() - 1], vertexes_[0]);
  if (curr.CrossesSegment(tmp)) {
    ++intersects;
  }
  for (uint32_t i = 0; i < vertexes_.size() - 1; ++i) {
    curr = Segment(vertexes_[i], vertexes_[i + 1]);
    if (curr.CrossesSegment(tmp)) {
      ++intersects;
    }
  }
  if (intersects % 2 == 0) {
    return !contains;
  }
  return contains;
}
bool Polygon::CrossesSegment(const Segment& seg) const {
  bool crosses = true;
  Segment tmp = Segment(vertexes_[vertexes_.size() - 1], vertexes_[0]);
  if (tmp.CrossesSegment(seg)) {
    return crosses;
  }
  for (uint32_t i = 0; i < vertexes_.size() - 1; ++i) {
    tmp = Segment(vertexes_[i], vertexes_[i + 1]);
    if (tmp.CrossesSegment(seg)) {
      return crosses;
    }
  }
  crosses = false;
  return crosses;
}

Circle::Circle(const Point& centre, int64_t radius) {
  centre_.x = centre.x;
  centre_.y = centre.y;
  radius_ = radius;
}

std::string Circle::ToString() {
  std::string circle_str =
      "Circle(" + centre_.ToString() + ", " + std::to_string(radius_) + ")";
  return circle_str;
}
IShape& Circle::Move(const Vector& vec) {
  centre_.x += vec.ReturnA().x;
  centre_.y += vec.ReturnA().y;
  return *this;
}
bool Circle::ContainsPoint(const Point& c) const {
  bool contains = true;
  int64_t dist = (centre_.x - c.x) * (centre_.x - c.x) +
                 (centre_.y - c.y) * (centre_.y - c.y);
  if (dist <= radius_ * radius_) {
    return contains;
  }
  contains = false;
  return contains;
}
bool Circle::CrossesSegment(const Segment& seg) const {
  bool crosses = true;
  double dist = seg.FindDistWithPoint(centre_);
  int64_t dist1 =
      (centre_.x - seg.ReturnBegin().x) * (centre_.x - seg.ReturnBegin().x) +
      (centre_.y - seg.ReturnBegin().y) * (centre_.y - seg.ReturnBegin().y);
  int64_t dist2 =
      (centre_.x - seg.ReturnEnd().x) * (centre_.x - seg.ReturnEnd().x) +
      (centre_.y - seg.ReturnEnd().y) * (centre_.y - seg.ReturnEnd().y);

  if (dist <= static_cast<double>(radius_) &&
      (dist1 >= radius_ * radius_ || dist2 >= radius_ * radius_)) {
    return crosses;
  }
  crosses = false;
  return crosses;
}

IShape* Circle::Clone() const {
  auto* cp = new Circle(centre_, radius_);
  return cp;
}

}  // namespace Geometry
template <class SmartPtrT>
void Delete(const SmartPtrT& ptr) {}

template <class T>
void Delete(T* ptr) {
  delete ptr;
}

void CheckFunctions(const Geometry::IShape* shape_ptr,
                    const Geometry::Point& point_a,
                    const Geometry::Point& point_b) {
  std::cout << "Given shape "
            << (shape_ptr->ContainsPoint(point_a) ? "contains"
                                                  : "does not contain")
            << " point A\n";

  const auto kSegmentAb = Geometry::Segment(point_a, point_b);
  std::cout << "Given shape "
            << (shape_ptr->CrossesSegment(kSegmentAb) ? "crosses"
                                                      : "does not cross")
            << " segment AB\n";

  const auto kVectorAb = point_b - point_a;
  const auto kClonedShapePtr =
      shape_ptr->Clone();  // may return either raw or smart pointer
  std::cout << kClonedShapePtr->Move(kVectorAb).ToString();

  Delete(kClonedShapePtr);  // raw pointer compatibility
}

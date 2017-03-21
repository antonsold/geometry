#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

const double EPS = 1E-7, PI = 3.14159265359; //машинный ноль

double det(double a, double b, double c, double d){
    return a * d - c * b;
}

class Point;

class Vector{
private:
    double x;
    double y;
public:
    Vector(): x(0), y(0){}
    Vector(double a, double b): x(a), y(b){}
    Vector(double x1, double y1, double x2, double y2): x(x2 - x1), y(y2 - y1){}
    Vector(const Vector& v):x(v.x), y(v.y){}
    Vector(const Point&, const Point&);
    double getX() const{
        return this->x;
    }
    double getY() const{
        return this->y;
    }
    double length() const{
        return sqrt(x * x + y * y);
    }
    Vector operator+(const Vector& v) const{
        return Vector(this->x + v.x, this->y + v.y);
    }
    Vector operator-(const Vector& v) const{
        return Vector(this->x - v.x, this->y - v.y);
    }
    Vector operator-() const{
        return Vector(-this->x, -this->y);
    }
    Vector operator*(const double k) const{
        return Vector(k * this->x, k * this->y);
    }
    Vector operator/(const double k) const{
        return Vector(this->x / k, this->y / k);
    }
    double operator*(const Vector& v) const {
        return this->x * v.x + this->y * v.y;
    }
    double vectorproduct(const Vector& v) const {
        return this->x * v.y - this->y * v.x;
    }
    friend ostream& operator<<(ostream &os, const Vector& v){
        os << v.x << ' ' << v.y;
        return os;
    }
};

class Segment;
class Line;
class Ray;
class Polygon;

class obj{
public:
    virtual void move(const Vector &) = 0;
    virtual bool contains(const Point &) const = 0;
    virtual bool intersects(const Segment &) const = 0;
};

class Point: public obj{
private:
    double x;
    double y;
public:
    Point(): x(0), y(0){}
    Point(double a, double b): x(a), y(b){}
    Point(const Point& p): x(p.x), y(p.y){}
    double getX() const{
        return this->x;
    }
    double getY() const{
        return this->y;
    }
    void move(const Vector& v){
        this->x += v.getX();
        this->y += v.getY();
    }
    bool contains(const Point &p) const{
        return abs(this->x - p.x) < EPS && abs(this->y - p.y) < EPS;
    }
    bool intersects(const Segment &s) const;
    double distance(const Point& p) const{
        return sqrt((p.x - this->x) * (p.x - this->x) + (p.y - this->y) * (p.y - this->y));
    }
    Point closestPoint(const Segment &s) const;
    friend ostream& operator<<(ostream& os, const Point& p){
        os << '(' << p.x << ", " << p.y << ')';
        return os;
    }
    friend istream& operator>>(istream& is, Point& p){
        is >> p.x >> p.y;
        return is;
    }
};

class Segment: public obj{
private:
    Point start;
    Point end;
public:
    Segment(double x1, double y1, double x2, double y2):start(Point(x1, y1)), end(Point(x2, y2)){};
    Segment(const Segment& s): start(s.start), end(s.end){}
    Segment(const Point& p1, const Point& p2):start(p1), end(p2){};
    const Point& getStart() const{
        return this->start;
    }
    const Point& getEnd() const{
        return this->end;
    }
    void move(const Vector& v){
        this->start.move(v);
        this->end.move(v);
    }
    bool contains(const Point& p) const{
        Vector AB(this->start, this->end);
        Vector CA(p, this->start);
        Vector CB(p, this->end);
        return abs(AB.vectorproduct(CA)) < EPS && CA * CB <= EPS;
    }
    bool intersects(const Segment& s) const{
        if (s.contains(this->start) || s.contains(this->end) || this->contains(s.getStart()) || this->contains(
                s.getEnd()))
            return true;
        else{
            Vector A1A2(this->start, this->end);
            Vector A1B1(this->start, s.start);
            Vector A1B2(this->start, s.end);
            return A1A2.vectorproduct(A1B1) * A1A2.vectorproduct(A1B2) < 0;
        }
    }
    double distanceToPoint(const Point& p) const;
};

bool Point::intersects(const Segment &s) const {
    Vector AB(s.getStart(), s.getEnd());
    Vector CA(*this, s.getStart());
    Vector CB(*this, s.getEnd());
    return abs(AB.vectorproduct(CA)) < EPS && CA * CB <= EPS;
}

Vector::Vector(const Point& p1, const Point &p2){
    this->x = p2.getX() - p1.getX();
    this->y = p2.getY() - p1.getY();
}

class Line: public obj{
private:
    double a;
    double b;
    double c;
public:
    Line(double _a, double _b, double _c): a(_a), b(_b), c(_c){}
    Line(const Point& p1, const Point& p2){
        this->a = p1.getY() - p2.getY();
        this->b = p2.getX() - p1.getX();
        this->c = det(p1.getX(), p1.getY(), p2.getX(), p2.getY());
    }
    Line(const Point& p, const Vector& v): a(-v.getY()), b(v.getX()), c(v.getY() * p.getX() - v.getX() * p.getY()){};
    Line(const Line& l): a(l.a), b(l.b), c(l.c){}
    Vector getVector(){
        Vector res(-this->b, this->a);
        return res;
    }
    double getA(){
        return this->a;
    }
    double getB(){
        return this->b;
    }
    double getC(){
        return this->c;
    }
    Point getPoint() const{ //получить точку на прямой
        if (a != 0)
            return Point(-c / a, 0);
        else
            return Point(0, -c / b);
    }
    void move(const Vector& v){
        this->c += this->a * v.getY() + this->b * v.getX();
    }
    bool contains(const Point& p) const{ //проверил на informatics
        return abs(this->a * p.getX() + this->b * p.getY() + this->c) < EPS;
    }
    bool intersects(const Segment& s) const{
        double x1, x2, y1, y2;
        if (this->a != 0){
            y1 = s.getStart().getY();
            y2 = s.getEnd().getY();
            x1 = (c - b * y1) / a;
            x2 = (c - b * y2) / a;
        }
        else{
            x1 = s.getStart().getX();
            x2 = s.getEnd().getX();
            y1 = (c - a * x1) / b;
            y2 = (c - a * x2) / b;
        }
        Segment pr(Point(x1, y1), Point(x2, y2));
        return s.intersects(pr);
    }
    bool intersects(const Line& line) const{
        double d = det(this->a, this->b, line.a, line.b);
        return abs(d) >= EPS;
    }
    double distanceToPoint(const Point &p) const{
        return abs(this->a * p.getX() + this->b * p.getY() + c) / sqrt(this->a * this->a + this->b * this->b);
    }
    double distanceToLine(const Line& line) const{ //только для параллельных или совпадающих
        Point p = line.getPoint();
        return this->distanceToPoint(p);
    }
    Point intersectPoint(const Line& line) const{ //только для пересекающихся
        double d = det(this->a, this->b, line.a, line.b);
        double x = -det(this->c, this->b, line.c, line.b) / d;
        double y = -det(this->a, this->c, line.a, line.c) / d;
        return Point(x, y);
    }
};

double Segment::distanceToPoint(const Point &p) const{ //проверил на informatics, работает
    Vector AP(this->start, p);
    Vector BP(this->end, p);
    Vector AB(this->start, this->end);
    if (AP * AB < EPS || BP * AB > -EPS)
        return min(p.distance(this->start), p.distance(this->end));
    else
        return (abs(AB.vectorproduct(AP))/ abs(AB.length()));
}

Point Point::closestPoint(const Segment &s) const{
    Vector AP(s.getStart(), *this);
    Vector BP(s.getEnd(), *this);
    Vector AB(s.getStart(), s.getEnd());
    if (AP * AB < 0)
        return s.getStart();
    else if (BP * AB > 0)
        return s.getEnd();
    else{
        double x, y;
        double xp = this->x, yp = this->y, x1 = s.getStart().x, y1 = s.getStart().y, x2 = s.getEnd().x,
                y2 = s.getEnd().y;
        if (y1 == y2){
            x = xp;
            y = y1;
        }
        else if (x1 == x2){
            x = x1;
            y = yp;
        }
        else{
            x = (x1 * (y2 - y1) * (y2 - y1) + xp * (x2 - x1) * (x2 - x1) + (x2 - x1) * (y2 - y1) * (yp - y1)) /
                    ((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)); //формула с википедии, вроде работает
            y = y1 + (x - x1) * (y2 - y1) / (x2 - x1);
        }
        return Point(x, y);
    }
}

class Ray: public obj{
private:
    Point start;
    Vector a;
public:
    Ray(const Point &p, const Vector &v):start(p), a(v){}
    Ray(const Point &p1, const Point &p2): start(p1), a(Vector(p1, p2)){}
    Ray(const Ray& r): start(r.start), a(r.a){}
    Point getPoint() const{
        Point p = this->start;
        p.move(this->a);
        return p;
    }
    void move(const Vector& v){
        this->start.move(v);
    }
    bool contains(const Point &p) const{ //проверил на informatics, работает
        Vector OP(this->start, p);
        return abs(OP.vectorproduct(this->a)) < EPS && OP * this->a > -EPS;
    }
    bool intersects(const Segment &s) const{
        Line rayLine(this->start, this->a);
        Line segmentLine(s.getStart(), s.getEnd());
        if (rayLine.intersects(segmentLine)){
            Point p = rayLine.intersectPoint(segmentLine);
            return s.contains(p) && this->contains(p);
        }
        else
            return false;
    }
    double distanceToPoint(const Point &p) const{
        Point B = this->getPoint();
        Vector OP(this->start, p);
        if (OP * this->a < EPS)
            return p.distance(this->start);
        else
            return abs(OP.vectorproduct(this->a) / OP.length());
    }
};

struct Vertex{
    Point p;
    Vertex* next;
};

class Polygon: public obj{
private:
    int size;
    Vertex *head, *tail;
public:
    Polygon(): size(0), head(NULL), tail(NULL){};
    ~Polygon(){
        while(size != 0){
            Vertex* temp = head->next;
            delete head;
            head = temp;
            size--;
        }
    }
    void add(Point p){
        size++;
        Vertex* temp = new Vertex;
        temp->p = p;
        temp->next = head;
        if (head != NULL){
            tail->next = temp;
            tail = temp;
        }
        else{
            tail = temp;
            head = temp;
        }
    }
    void move(const Vector &v){
        Vertex* tmp = head;
        do{
            tmp->p.move(v);
            tmp = tmp->next;
        }
        while (tmp != tail);
    }
    void print(){
        Vertex* tmp = head;
        for (int i = 0; i < this->size; i++){
            cout << tmp->p;
            tmp = tmp->next;
        }
    }
    bool contains(const Point &p) const;
    bool intersects(const Segment &s) const;
};

bool Polygon::contains(const Point &p) const { //кольцевой связный список
    Vertex* tmp = head;
    Point lastP = this->tail->p;
    do{ //проверка, лежит ли точка на границе
        Segment s(lastP, tmp->p);
        lastP = tmp->p;
        if (s.contains(p))
            return true;
        tmp = tmp->next;
    }
    while(tmp != head);
    lastP = head->p;
    double angleSum = 0;
    do { //метод суммирования углов
        double part = 0;
        tmp = tmp->next;
        Vector vLast(p, lastP);
        Vector vCur(p, tmp->p);
        part = acos(vLast * vCur / (vLast.length() * vCur.length()));
        cout << part << endl;
        if (det(vLast.getX(), vLast.getY(), vCur.getX(), vCur.getY()) > -EPS)
            angleSum += part;
        else
            angleSum -= part;
        lastP = tmp->p;
    }
    while(tmp != head);
    return (abs(angleSum /(2 * PI) - 1) < EPS || abs(angleSum /(2 * PI) + 1) < EPS);
}

bool Polygon::intersects(const Segment &s) const {
    Vertex* tmp = head;
    Point lastP = tmp->p;
    bool flag = false;
    do{
        tmp = tmp->next;
        Segment edge(lastP, tmp->p);
        if (s.intersects(edge))
            flag = true;
        lastP = tmp->p;
    }
    while(!flag && tmp != head);
    return flag;
}

int main() {
    //решение задачи о принадлежности точки многоугольнику
    int n;
    Point cur;
    cin >> n;
    Point pt;
    cin >> pt;
    Polygon p;
    for (int i = 0; i < n; i++) {
        cin >> cur;
        p.add(cur);
    }
    if (p.contains(pt))
        cout << "YES";
    else
        cout << "NO";
    return 0;
}
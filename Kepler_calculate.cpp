#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>

using namespace std;
using namespace Eigen;

// 定义π常量
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// 常量定义
const double TOLERANCE = 1e-8;
const int MAX_ITERATIONS = 100;
const double EARTH_MU = 398600.4418;  // 地球引力常数 (km^3/s^2)

// 轨道根数结构体
struct OrbitalElements {
    double a;        // 半长轴 (km)
    double e;        // 偏心率
    double i;        // 轨道倾角 (rad)
    double omega;    // 近地点幅角 (rad)
    double Omega;    // 升交点赤经 (rad)
    double f;        // 真近点角 (rad)

    OrbitalElements() : a(0), e(0), i(0), omega(0), Omega(0), f(0) {} //初始化结构体

    // 规范化所有角度真近点角、升交点赤经、近地点幅角到 [0, 2π) 范围
    void normalizeAngles() {
        omega = fmod(omega, 2 * M_PI);
        if (omega < 0) omega += 2 * M_PI;

        Omega = fmod(Omega, 2 * M_PI);
        if (Omega < 0) Omega += 2 * M_PI;

        f = fmod(f, 2 * M_PI);
        if (f < 0) f += 2 * M_PI;
    }
};

// 位置速度状态向量结构体
struct StateVector {
    Vector3d position;  // 位置矢量 (km)
    Vector3d velocity;  // 速度矢量 (km/s)

    StateVector() : position(Vector3d::Zero()), velocity(Vector3d::Zero()) {}
    StateVector(const Vector3d& pos, const Vector3d& vel) : position(pos), velocity(vel) {}
};

// 旋转矩阵
namespace Rotation {
    Matrix3d Rx(double angle) {
        Matrix3d mat;
        mat << 1, 0, 0,
            0, cos(angle), -sin(angle),
            0, sin(angle), cos(angle);
        return mat;
    }

    Matrix3d Rz(double angle) {
        Matrix3d mat;
        mat << cos(angle), -sin(angle), 0,
            sin(angle), cos(angle), 0,
            0, 0, 1;
        return mat;
    }
}

// 开普勒方程求解器 用于迭代求解偏近点角
class KeplerSolver {
public:
    static double solve(double M, double e) {
        double E = M;  // 初始值
        int iterations = 0;
        double delta;

        do {
            double f = E - e * sin(E) - M;
            double f_prime = 1 - e * cos(E);
            delta = f / f_prime;
            E -= delta;
            iterations++;
        } while (fabs(delta) > TOLERANCE && iterations < MAX_ITERATIONS);

        if (iterations >= MAX_ITERATIONS) {
            throw runtime_error("Kepler equation did not converge");
        }

        return E;
    }
};

// Kepler->PV && PV->Kepler 轨道转换类 用于实现 由状态向量到轨道根数 和 由轨道根数到状态向量 的计算
class OrbitConverter {
private:
    double mu_;
public:
    explicit OrbitConverter(double mu = EARTH_MU) : mu_(mu) {}  //使用默认参数的构造函数 

    // 获取引力常数
    double getMu() const { return mu_; }

    // 从位置、速度计算轨道根数
    OrbitalElements stateToElements(const StateVector& state) const {
        OrbitalElements elements;

        Vector3d r = state.position;
        Vector3d v = state.velocity;

        // 计算角动量矢量
        Vector3d h = r.cross(v);
        double h_norm = h.norm();

        // 计算轨道倾角
        elements.i = acos(h(2) / h_norm);

        // 计算升交点赤经
        Vector3d n(-h(1), h(0), 0);
        double n_norm = n.norm();
        if (n_norm < TOLERANCE) {
            elements.Omega = 0; // 赤道轨道
        }
        else {
            elements.Omega = atan2(n(1), n(0));
            if (elements.Omega < 0) elements.Omega += 2 * M_PI;
        }

        // 计算偏心率矢量
        Vector3d e_vec = v.cross(h) / mu_ - r.normalized();
        elements.e = e_vec.norm();

        // 计算近地点幅角
        if (n_norm > TOLERANCE) { // 非赤道轨道
            Vector3d node_vec = n.normalized();
            Vector3d perp_vec = h.normalized().cross(node_vec);

            double e_dot_node = e_vec.dot(node_vec);
            double e_dot_perp = e_vec.dot(perp_vec);
            elements.omega = atan2(e_dot_perp, e_dot_node);
        }
        else { // 赤道轨道
            elements.omega = atan2(e_vec(1), e_vec(0));
        }
        if (elements.omega < 0) elements.omega += 2 * M_PI;

        // 计算真近点角 (使用atan2确保正确象限)
        double p = (h_norm * h_norm) / mu_;
        double r_norm = r.norm();

        if (elements.e > TOLERANCE) { // 椭圆轨道
            double cos_f = (p / r_norm - 1) / elements.e;
            double sin_f = (r.dot(v)) / r_norm * (p / (elements.e * h_norm));
            elements.f = atan2(sin_f, cos_f);
        }
        else { // 圆轨道
            Vector3d node_vec(cos(elements.Omega), sin(elements.Omega), 0);
            Vector3d perp_vec = h.normalized().cross(node_vec);
            elements.f = atan2(r.dot(perp_vec), r.dot(node_vec));
        }
        if (elements.f < 0) elements.f += 2 * M_PI;

        // 计算半长轴
        double energy = v.squaredNorm() / 2 - mu_ / r_norm;
        elements.a = -mu_ / (2 * energy);

        elements.normalizeAngles();
        return elements;
    }

    // 从轨道根数计算状态向量
    StateVector elementsToState(const OrbitalElements& elements) const {
        // 计算轨道坐标系中的位置和速度
        double p = elements.a * (1 - elements.e * elements.e);
        double r = p / (1 + elements.e * cos(elements.f));

        Vector3d r_orb(r * cos(elements.f), r * sin(elements.f), 0);
        double sqrt_mu_p = sqrt(mu_ / p);
        Vector3d v_orb(-sqrt_mu_p * sin(elements.f),
            sqrt_mu_p * (elements.e + cos(elements.f)),
            0);

        // 构建旋转矩阵：轨道坐标系 -> 惯性坐标系
        Matrix3d R = Rotation::Rz(elements.Omega) *
            Rotation::Rx(elements.i) *
            Rotation::Rz(elements.omega);

        // 转换到惯性系
        Vector3d r_inertial = R * r_orb;
        Vector3d v_inertial = R * v_orb;

        return StateVector(r_inertial, v_inertial);
    }
};

// 轨道预测
class OrbitPredictor {
private:
    OrbitConverter converter_; //私有变量为轨道转换类

public:
    explicit OrbitPredictor(double mu = EARTH_MU) : converter_(mu) {} // 默认参数地球引力常数构造函数

    // 获取引力常数
    double getMu() const { return converter_.getMu(); }

    // 预报未来时刻的状态 参数为初始状态向量和间隔时间deltaTime
    StateVector predict(const StateVector& initialState, double deltaTime) const {
        // 转换为轨道根数
        OrbitalElements elements = converter_.stateToElements(initialState);

        // 计算初始偏近点角
        double E0 = 2 * atan(sqrt((1 - elements.e) / (1 + elements.e)) * tan(elements.f / 2));
        E0 = fmod(E0, 2 * M_PI);
        if (E0 < 0) E0 += 2 * M_PI;

        // 计算初始平近点角
        double M0 = E0 - elements.e * sin(E0);

        // 计算平均角速度
        double n = sqrt(converter_.getMu() / pow(elements.a, 3));

        // 计算未来平近点角并规范化
        double Mt = M0 + n * deltaTime;
        Mt = fmod(Mt, 2 * M_PI);
        if (Mt < 0) Mt += 2 * M_PI;

        // 解算未来偏近点角 开普勒方程迭代求解 并规范化角度
        double Et = KeplerSolver::solve(Mt, elements.e);
        Et = fmod(Et, 2 * M_PI);
        if (Et < 0) Et += 2 * M_PI;

        // 计算未来真近点角
        elements.f = 2 * atan2(sqrt(1 + elements.e) * sin(Et / 2),
            sqrt(1 - elements.e) * cos(Et / 2));
        if (elements.f < 0) elements.f += 2 * M_PI;

        // 转换回状态向量
        return converter_.elementsToState(elements);
    }
};
int main() {
    try {
        // 创建初始位置、速度 
        Vector3d position(7076.152744, -487.8211097, -17.539112);  // 本人学号2023302141161 
        Vector3d velocity(-0.053864, -1.041859, 7.461029);
        StateVector initialState(position, velocity);// 初始化状态向量

        // 位置、速度预测
        OrbitPredictor predictor;

        // 预报1小时后的状态
        double deltaTime = 3600.0;  // 1小时
        StateVector futureState = predictor.predict(initialState, deltaTime);

        // 输出结果
        cout << "Initial position: " << initialState.position.transpose() << " km" << endl;
        cout << "Initial velocity: " << initialState.velocity.transpose() << " km/s" << endl;
        cout << "Future position: " << futureState.position.transpose() << " km" << endl;
        cout << "Future velocity: " << futureState.velocity.transpose() << " km/s" << endl;

        // 验证轨道根数转换的准确性
        OrbitConverter converter;
        OrbitalElements elements = converter.stateToElements(initialState);
        StateVector reconstructedState = converter.elementsToState(elements);

        // 用实际位置和预测位置的直线距离代表预测误差
        cout << "\nReconstruction error: "
            << (initialState.position - reconstructedState.position).norm() << " km" << endl;  

    }
    catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }

    return 0;
}
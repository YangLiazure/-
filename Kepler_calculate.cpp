#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>

using namespace std;
using namespace Eigen;

// ����г���
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ��������
const double TOLERANCE = 1e-8;
const int MAX_ITERATIONS = 100;
const double EARTH_MU = 398600.4418;  // ������������ (km^3/s^2)

// ��������ṹ��
struct OrbitalElements {
    double a;        // �볤�� (km)
    double e;        // ƫ����
    double i;        // ������ (rad)
    double omega;    // ���ص���� (rad)
    double Omega;    // ������ྭ (rad)
    double f;        // ������ (rad)

    OrbitalElements() : a(0), e(0), i(0), omega(0), Omega(0), f(0) {} //��ʼ���ṹ��

    // �淶�����нǶ������ǡ�������ྭ�����ص���ǵ� [0, 2��) ��Χ
    void normalizeAngles() {
        omega = fmod(omega, 2 * M_PI);
        if (omega < 0) omega += 2 * M_PI;

        Omega = fmod(Omega, 2 * M_PI);
        if (Omega < 0) Omega += 2 * M_PI;

        f = fmod(f, 2 * M_PI);
        if (f < 0) f += 2 * M_PI;
    }
};

// λ���ٶ�״̬�����ṹ��
struct StateVector {
    Vector3d position;  // λ��ʸ�� (km)
    Vector3d velocity;  // �ٶ�ʸ�� (km/s)

    StateVector() : position(Vector3d::Zero()), velocity(Vector3d::Zero()) {}
    StateVector(const Vector3d& pos, const Vector3d& vel) : position(pos), velocity(vel) {}
};

// ��ת����
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

// �����շ�������� ���ڵ������ƫ�����
class KeplerSolver {
public:
    static double solve(double M, double e) {
        double E = M;  // ��ʼֵ
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

// Kepler->PV && PV->Kepler ���ת���� ����ʵ�� ��״̬������������� �� �ɹ��������״̬���� �ļ���
class OrbitConverter {
private:
    double mu_;
public:
    explicit OrbitConverter(double mu = EARTH_MU) : mu_(mu) {}  //ʹ��Ĭ�ϲ����Ĺ��캯�� 

    // ��ȡ��������
    double getMu() const { return mu_; }

    // ��λ�á��ٶȼ���������
    OrbitalElements stateToElements(const StateVector& state) const {
        OrbitalElements elements;

        Vector3d r = state.position;
        Vector3d v = state.velocity;

        // ����Ƕ���ʸ��
        Vector3d h = r.cross(v);
        double h_norm = h.norm();

        // ���������
        elements.i = acos(h(2) / h_norm);

        // ����������ྭ
        Vector3d n(-h(1), h(0), 0);
        double n_norm = n.norm();
        if (n_norm < TOLERANCE) {
            elements.Omega = 0; // ������
        }
        else {
            elements.Omega = atan2(n(1), n(0));
            if (elements.Omega < 0) elements.Omega += 2 * M_PI;
        }

        // ����ƫ����ʸ��
        Vector3d e_vec = v.cross(h) / mu_ - r.normalized();
        elements.e = e_vec.norm();

        // ������ص����
        if (n_norm > TOLERANCE) { // �ǳ�����
            Vector3d node_vec = n.normalized();
            Vector3d perp_vec = h.normalized().cross(node_vec);

            double e_dot_node = e_vec.dot(node_vec);
            double e_dot_perp = e_vec.dot(perp_vec);
            elements.omega = atan2(e_dot_perp, e_dot_node);
        }
        else { // ������
            elements.omega = atan2(e_vec(1), e_vec(0));
        }
        if (elements.omega < 0) elements.omega += 2 * M_PI;

        // ���������� (ʹ��atan2ȷ����ȷ����)
        double p = (h_norm * h_norm) / mu_;
        double r_norm = r.norm();

        if (elements.e > TOLERANCE) { // ��Բ���
            double cos_f = (p / r_norm - 1) / elements.e;
            double sin_f = (r.dot(v)) / r_norm * (p / (elements.e * h_norm));
            elements.f = atan2(sin_f, cos_f);
        }
        else { // Բ���
            Vector3d node_vec(cos(elements.Omega), sin(elements.Omega), 0);
            Vector3d perp_vec = h.normalized().cross(node_vec);
            elements.f = atan2(r.dot(perp_vec), r.dot(node_vec));
        }
        if (elements.f < 0) elements.f += 2 * M_PI;

        // ����볤��
        double energy = v.squaredNorm() / 2 - mu_ / r_norm;
        elements.a = -mu_ / (2 * energy);

        elements.normalizeAngles();
        return elements;
    }

    // �ӹ����������״̬����
    StateVector elementsToState(const OrbitalElements& elements) const {
        // ����������ϵ�е�λ�ú��ٶ�
        double p = elements.a * (1 - elements.e * elements.e);
        double r = p / (1 + elements.e * cos(elements.f));

        Vector3d r_orb(r * cos(elements.f), r * sin(elements.f), 0);
        double sqrt_mu_p = sqrt(mu_ / p);
        Vector3d v_orb(-sqrt_mu_p * sin(elements.f),
            sqrt_mu_p * (elements.e + cos(elements.f)),
            0);

        // ������ת���󣺹������ϵ -> ��������ϵ
        Matrix3d R = Rotation::Rz(elements.Omega) *
            Rotation::Rx(elements.i) *
            Rotation::Rz(elements.omega);

        // ת��������ϵ
        Vector3d r_inertial = R * r_orb;
        Vector3d v_inertial = R * v_orb;

        return StateVector(r_inertial, v_inertial);
    }
};

// ���Ԥ��
class OrbitPredictor {
private:
    OrbitConverter converter_; //˽�б���Ϊ���ת����

public:
    explicit OrbitPredictor(double mu = EARTH_MU) : converter_(mu) {} // Ĭ�ϲ������������������캯��

    // ��ȡ��������
    double getMu() const { return converter_.getMu(); }

    // Ԥ��δ��ʱ�̵�״̬ ����Ϊ��ʼ״̬�����ͼ��ʱ��deltaTime
    StateVector predict(const StateVector& initialState, double deltaTime) const {
        // ת��Ϊ�������
        OrbitalElements elements = converter_.stateToElements(initialState);

        // �����ʼƫ�����
        double E0 = 2 * atan(sqrt((1 - elements.e) / (1 + elements.e)) * tan(elements.f / 2));
        E0 = fmod(E0, 2 * M_PI);
        if (E0 < 0) E0 += 2 * M_PI;

        // �����ʼƽ�����
        double M0 = E0 - elements.e * sin(E0);

        // ����ƽ�����ٶ�
        double n = sqrt(converter_.getMu() / pow(elements.a, 3));

        // ����δ��ƽ����ǲ��淶��
        double Mt = M0 + n * deltaTime;
        Mt = fmod(Mt, 2 * M_PI);
        if (Mt < 0) Mt += 2 * M_PI;

        // ����δ��ƫ����� �����շ��̵������ ���淶���Ƕ�
        double Et = KeplerSolver::solve(Mt, elements.e);
        Et = fmod(Et, 2 * M_PI);
        if (Et < 0) Et += 2 * M_PI;

        // ����δ��������
        elements.f = 2 * atan2(sqrt(1 + elements.e) * sin(Et / 2),
            sqrt(1 - elements.e) * cos(Et / 2));
        if (elements.f < 0) elements.f += 2 * M_PI;

        // ת����״̬����
        return converter_.elementsToState(elements);
    }
};
int main() {
    try {
        // ������ʼλ�á��ٶ� 
        Vector3d position(7076.152744, -487.8211097, -17.539112);  // ����ѧ��2023302141161 
        Vector3d velocity(-0.053864, -1.041859, 7.461029);
        StateVector initialState(position, velocity);// ��ʼ��״̬����

        // λ�á��ٶ�Ԥ��
        OrbitPredictor predictor;

        // Ԥ��1Сʱ���״̬
        double deltaTime = 3600.0;  // 1Сʱ
        StateVector futureState = predictor.predict(initialState, deltaTime);

        // ������
        cout << "Initial position: " << initialState.position.transpose() << " km" << endl;
        cout << "Initial velocity: " << initialState.velocity.transpose() << " km/s" << endl;
        cout << "Future position: " << futureState.position.transpose() << " km" << endl;
        cout << "Future velocity: " << futureState.velocity.transpose() << " km/s" << endl;

        // ��֤�������ת����׼ȷ��
        OrbitConverter converter;
        OrbitalElements elements = converter.stateToElements(initialState);
        StateVector reconstructedState = converter.elementsToState(elements);

        // ��ʵ��λ�ú�Ԥ��λ�õ�ֱ�߾������Ԥ�����
        cout << "\nReconstruction error: "
            << (initialState.position - reconstructedState.position).norm() << " km" << endl;  

    }
    catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }

    return 0;
}
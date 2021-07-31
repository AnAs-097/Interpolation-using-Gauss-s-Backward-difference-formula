// NA assignment 3.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <conio.h>
using namespace std;

class NewBackDiff
{
public:
    NewBackDiff();
    void input();
    void result();
    void get_1st_der();
    void get_2nd_der();
private:
    int degree, values, actual_degree;
    float sigma[10], xp, p, x[10], fx[10], h;
};

NewBackDiff::NewBackDiff()
{
    cout << "\n\t---------GAUSS'S BACKWARD DIFFERENCE FORMULA--------\n\n";
    cout << "\n\tfp = fo + p.&f(-1/2) + 1/2(p^2 +p)&^2fo + 1/6(p^3-p)&^3f(-1/2) + 1/24(p^4 - 2p^3 - p^2 + 2p)&^4fo + ....\n\n";

    degree = values = xp = actual_degree = 0;
    p = -2;
    for (int i = 0; i < 10; i++)
        sigma[i] = x[i] = fx[i] = 0.0;
}

void NewBackDiff::input()
{
    cout << "How many values of X?\t";
    cin >> values;
    cout << "Upto what power of Sigma?\t";
    cin >> degree;
    cout << "\n Value of Xp:\t";
    cin >> xp;
    for (int i = 0; i != values; i++)
    {
        cout << "\nEnter X" << i + 1 << ":\t";
        cin >> x[i]; cout << "Enter F(" << i + 1 << "):\t";
        cin >> fx[i];
    }
}

void NewBackDiff::result()
{
    int temp = -1, origin = 0;
    cout << "\nX\t";
    for (int i = 0; i != values; i++)
        cout << "\t" << x[i];
    cout << endl;
    cout << "\nF(x)\t";
    for (int i = 0; i != values; i++)
        cout << "\t" << fx[i];

    cout << endl;
    h = x[1] - x[0];
    cout << endl << "h: " << h << endl;
    for (int j = values - 1; j >= 0 && (p < 0 || p>1); j--)
    {
        p = (xp - x[j]) / h;
        temp = j;
    }
    cout << "\nValue of P is :\t" << p << "\n";
    origin = temp;
    actual_degree = 1;

    //Printing Difference table
    for (int j = values; actual_degree <= degree && j > 1; actual_degree++, j--)
    {
        cout << "\n\nsigma power" << actual_degree << ":\t";
        for (int k = 0; k < j - 1; k++)
        {
            fx[k] = fx[k + 1] - fx[k];
            cout << fx[k] << "\t";
        }
        if (actual_degree % 2 == 1) {

            sigma[actual_degree - 1] = fx[temp - 1];

        }
        else if (actual_degree % 2 == 0) {
            sigma[actual_degree - 1] = fx[temp - 1];
            temp = temp - 1;
        }

    }
    get_1st_der();
    get_2nd_der();
}
//calcualte 1st derivative
/*void NewBackDiff::get_1st_der()
{
    float parray[] = { 1,2 * p + 1, 3 * p * p - 1,- (2 * p * p * p) - (3 * p * p) - p + 1 },
        div[] = { 1,2,6,12 },
        ans = 0;

    for (int i = 0; i < actual_degree; i++)
    {
        ans += sigma[i] * parray[i] / div[i];
    }
    ans = ans / h;
    cout << "\n\n\nf'(" << xp << "):\t" << ans;
}*/
void NewBackDiff::get_1st_der()
{
    cout <<endl<< "The formula for 1st order derivative: 1/h{&f(-1/2) + 1/2(2p + 1)&^2fo + 1/6(3p^2 - 1)&^3f(-1/2) + 1/12(2p^3 - 3p^2 - p + 1)&^4fo}"<<endl;

    float parray[] = { 1,2 * p + 1, 3 * p * p - 1, (2 * p * p * p) - (3 * p * p) - p + 1 },
        div[] = { 1,2,6,12 },
        ans = 0;

    for (int i = 0; i < actual_degree; i++)
    {
        ans += sigma[i] * parray[i] / div[i];
    }
    ans = ans / h;
    cout << "\n\n\nf'(" << xp << "):\t" << ans;
}
//calcualte 2nd derivative
void NewBackDiff::get_2nd_der()
{
    cout <<endl<<endl<< "The formula for 2nd order derivative: 1/h^2{&^2fo + p.&^3f(-1/2) + 1/12(6p^2 - 6p - 1)&^4fo}"<<endl;
    float parray[] = { 1 , p , (6 * p * p - 6 * p - 1) },
        div[] = { 1,1,12 },
        ans = 0;

    for (int i = 0; i < actual_degree; i++)
    {
        ans += sigma[i + 1] * parray[i] / div[i];
    }
    ans = ans / (h*h);
    cout << "\n\n\nf''(" << xp << "):\t" << ans;
}

int main()
{
    NewBackDiff obj;
    obj.input();
    obj.result();

    _getch();
}


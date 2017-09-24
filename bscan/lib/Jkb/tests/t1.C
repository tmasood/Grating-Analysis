#include <iostream.h>
#include <Jkb.h>

int main()
{
  int N;
  double x;

  cout << "Enter N:\n";
  cin >> N;

  cout << "Enter x:\n";
  cin >> x;

  cout << "Cheb ratio = " << dunounp1(N, x) << '\n' << endl;

}

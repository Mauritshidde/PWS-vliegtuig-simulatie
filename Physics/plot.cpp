#include <cmath>
#include <matplot/matplot.h>

int main()
{
      double x, y = 0;
      auto [x, y] = matplot::meshgrid(matplot::iota(0.0, 0.2, 2.0), matplot::iota(0.0, 0.2, 2.0));
      matplot::vector_2d u =
          matplot::transform(x, y, [](double x, double y)
                             { return cos(x) * y; });
      matplot::vector_2d v =
          matplot::transform(x, y, [](double x, double y)
                             { return sin(x) * y; });

      matplot::quiver(x, y, u, v);

      matplot::show();
      return 0;
}
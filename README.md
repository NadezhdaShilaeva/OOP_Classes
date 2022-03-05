# OOP_Classes
###### Лабораторная работа по программированию №1. 2 семестр.

## 1. Классы, наследование
Программа содержится в файле ___Shapes.cpp___.

В программе спроектированы и реализованы следующие классы:
1. ___Точка___
2. ___Ломаная___
3. ___Замкнутая ломаная___ (унаследована от ломаной)
4. ___Многоугольник___
5. ___Треугольник___ (унаследован от многоугольника)
6. ___Трапеция___ (унаследован от многоугольника)
7. ___Правильный многоугольник___ (унаследован от многоугольника)

Для каждого из классов реализованы следующие методы:
1. ___Конструкторы___
2. ___Конструктор копирования___
3. ___Оператор присваивания___
4. ___Расчет длины___ (для ломаной и замкнутой ломаной)
5. ___Расчет периметра___ (для многоугольника, треугольника, трапеции, правильного многоугольника)
6. ___Расчет площади___ (для многоугольника, треугольника, трапеции, правильного многоугольника)
7. ___Название фигуры___ (для многоугольника, треугольника, трапеции, правильного многоугольника)
8. ___Методы для доступа к значениям полей соответствующих классов___ (для точки, ломаной, многоугольника)

В функции ___main()___ продемонстрирован динамический полиморфизм. Здесь создается массив, элементами которого могут быть объекты таких калассов, как _многоугольник_, _треугольник_, _трапеция_, _правильный многоугольник_ (данные объекты вводятся пользователем), и далее выводятся соответствующие _название_, _периметр_ и _площадь_ для каждого элемента.

#### Формат входных данных:
Первым вводится натуральное число ___n___ - размер массива.
Далее вводятся n элементов массива, кажлый из которых описывается числом ___num_points___ - количество точек фигуры и соответствующими точками, заданными в координатах ___x___ и ___y___.

---
#### Примеры использования программы:
___Входные аргументы:___  
4  
4  
2 3 4 1 -2 -1 -1 2  
3  
2 1 7 8 4 5  
4  
0 3 3 0 0 -3 -3 0  
5  
3 4 5 11 12 8 9 5 5 6  
___Выходные данные:___  
Trapezoid 15.4775 12  
Triangle 17.3171 3  
Regular polygon 16.9706 18  
Polygon 26.0901 30  

## 2. Классы, перегрузка операторов
Программа содержится в файле ___Polynomial.cpp___.

В программе спроектирован и реализован класс для описания сущности многочлен (полином), раздела математики - Алгебра.

Для данного класса реализованы ___конструкторы___, ___конструктор копирования___, ___деструктор___, а также следующие операторы:
1. __=__
2. __==__, __!=__
3. __+__, __-__ (унарный и бинарный), __+=__, __-=__
4. __*__ (на многочлен и на число), __/__ (на число), __*=__, __/=__
5. __<<__, __>>__
6. __[]__ (для получения коэффициента i-го члена)

В функции ___main()___ на вход принимается два многочлена, программа выводит различные дествия, выполненые с ними.

#### Формат входных данных:
Вводятся два полинома - ___poly_left___ и ___poly_right___. Все коэффициенты необходимо прописать явно.

---
#### Примеры использования программы:
___Входные аргументы:___  
poly_left: 1x^3 + 1.5x^1  
poly_right: 0.5x^5 - 2.5x^4 + 3x^3 - 6.25x^0   
___Выходные данные:___  
poly_left != poly_right  
poly_left + poly_right = 0.5x^5 - 2.5x^4 + 4x^3 + 1.5x^1 - 6.25x^0    
poly_left - poly_right = -0.5x^5 + 2.5x^4 - 2x^3 + 1.5x^1 + 6.25x^0  
-poly_right = -0.5x^5 + 2.5x^4 - 3x^3 + 6.25x^0  
poly_left -= poly_right = -0.5x^5 + 2.5x^4 - 2x^3 + 1.5x^1 + 6.25x^0  
poly_left += poly_right = 1x^3 + 1.5x^1  
poly_left /= 2 = 0.5x^3 + 0.75x^1  
poly_left *= 2 = 1x^3 + 1.5x^1  
poly_left * poly_right = 0.5x^8 - 2.5x^7 + 3.75x^6 - 3.75x^5 + 4.5x^4 - 6.25x^3 - 9.375x^1  
poly_left degrees: 0 -9.375 0 -6.25 4.5 -3.75 3.75 -2.5 0.5  
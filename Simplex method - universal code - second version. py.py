# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 19:35:10 2025

@author: Stoqn
"""

#from Rationals import Rational
from sympy import symbols, Rational, simplify

M_VALUE = 10000
M = symbols("M")

def read_coefficients(num, text):
    print(f"Въведете коефициентите за {text}:")
    return [Rational(input(f"Коефициент {i + 1}: ")) for i in range(num)]

def read_constraints(num_vars, num_constraints):
    constraints = []
    for i in range(num_constraints):
        print(f"Въведете коефициентите за ограничение {i + 1}:")
        coeffs = [Rational(input(f"Коефициент {j + 1}: ")) for j in range(num_vars)]
        sign = input("Въведете знака на ограничението (<=, =, >=): ")
        rhs = Rational(input("Въведете дясната страна на ограничението: "))
        while rhs < 0:
            print('\n', 'Дясната страна на ограниченията трябва да бъде неотрицателна. Моля умножете ограничението по -1.', '\n')
            print(f"Въведете коефициентите за ограничение {i + 1}:")
            coeffs = [Rational(input(f"Коефициент {j + 1}: ")) for j in range(num_vars)]
            sign = input("Въведете знака на ограничението (<=, =, >=): ")
            rhs = Rational(input("Въведете дясната страна на ограничението: "))
        constraints.append((coeffs, sign, rhs))
    return constraints

def display_function(coeffs, objective, free_coeff):
    terms = [f"{coeff}*x{i + 1}" for i, coeff in enumerate(coeffs)]
    objective_str = "Максимизация" if objective == "max" else "Минимизация"
    
    # Присъединяване на свободния член само ако е различен от 0
    if free_coeff != 0:
        free_coeff_str = f" + {free_coeff}" if free_coeff > 0 else f" - {-free_coeff}"
    else:
        free_coeff_str = ""
    
    return f"{objective_str} на: " + " + ".join(terms).replace('+ -', '- ') + free_coeff_str 

def display_constraints(constraints):
    result = []
    for coeffs, sign, rhs in constraints:
        terms = [f"{coeff}*x{i + 1}" for i, coeff in enumerate(coeffs)]
        constraint = " + ".join(terms).replace('+ -', '- ') + f" {sign} {rhs}"
        result.append(constraint)
    return "\n".join(result)

def find_basis(constraints):
    """Намира базисните променливи в ограниченията."""
    num_constraints = len(constraints)
    num_vars = len(constraints[0][0]) if constraints else 0
    basis_vars = [-1] * num_constraints  # Индексите на базисните променливи

    for j in range(num_vars):
        count = 0
        row_index = -1
        for i in range(num_constraints):
            if constraints[i][0][j] == 1:
                count += 1
                row_index = i
            elif constraints[i][0][j] != 0:
                count = -1  # Ако има друга ненулева стойност, не е базис
                break
        if count == 1:
            basis_vars[row_index] = j

    return basis_vars

def canonicalize(coefficients, constraints, objective, free_coefficient):
    num_vars = len(coefficients)
    surplus_vars = 0
    
    has_equality_constraints = any(sign == "=" for _, sign, _ in constraints)
    has_inequality_constraints = any(sign == "<=" or sign == ">=" for _, sign, _ in constraints)
    count_equality_constraints = sum(1 for _, sign, _ in constraints if sign == "=")
    
    if has_equality_constraints and has_inequality_constraints:
        # Добавяне на излишъци променливи където е необходимо
        for i, (coeffs, sign, rhs) in enumerate(constraints):
            if sign == "<=":
                coeffs.extend([Rational(1) if j == surplus_vars else Rational(0) for j in range(len(constraints) - count_equality_constraints)])
                surplus_vars += 1
                coefficients.append(Rational(0))
            elif sign == ">=":
                coeffs.extend([Rational(-1) if j == surplus_vars else Rational(0) for j in range(len(constraints) - count_equality_constraints)])
                surplus_vars += 1
                coefficients.append(Rational(0))
            elif sign == "=":
                coeffs.extend([Rational(0)] * (len(constraints) - count_equality_constraints))
            constraints[i] = (coeffs, "=", rhs)
    
        # Актуализация на броя на променливите
        num_vars += surplus_vars
    
        # Проверка за базис
        basis_vars = find_basis(constraints)
    
        # Добавяне на изкуствени променливи, ако няма достатъчно базисни
        artificial_vars = []
        for i in range(len(constraints)):
            if basis_vars[i] == -1:  # Ако няма базисна променлива в това ограничение
                artificial_var_index = num_vars + len(artificial_vars)
                for j in range(len(constraints)):
                    constraints[j][0].append(Rational(1) if j == i else Rational(0))  # Добавяне на изкуствена променлива
                artificial_vars.append(artificial_var_index)
                if objective == "min":
                    coefficients.append(Rational(M_VALUE))
                elif objective == "max":
                    coefficients.append(Rational(-M_VALUE)) # Изкуствените променливи не влияят на целевата функция
    
        num_vars += len(artificial_vars)
    
        # Подготовка и извеждане на каноничния вид
        print("Каноничен вид на задачата:")
        print("Целева функция:")
        print_function(coefficients, num_vars, objective, free_coefficient)
    
        print("Ограничения:")
        for coeffs, sign, rhs in constraints:
            print_constraint(coeffs, rhs)
    else:    
        # Добавяне на излишъци променливи където е необходимо
        for i, (coeffs, sign, rhs) in enumerate(constraints):
            if sign == "<=":
                coeffs.extend([Rational(1) if j == surplus_vars else Rational(0) for j in range(len(constraints))])
                surplus_vars += 1
                coefficients.append(Rational(0))
            elif sign == ">=":
                coeffs.extend([Rational(-1) if j == surplus_vars else Rational(0) for j in range(len(constraints))])
                surplus_vars += 1
                coefficients.append(Rational(0))
            constraints[i] = (coeffs, "=", rhs)
    
        # Актуализация на броя на променливите
        num_vars += surplus_vars
    
        # Проверка за базис
        basis_vars = find_basis(constraints)
    
        # Добавяне на изкуствени променливи, ако няма достатъчно базисни
        artificial_vars = []
        for i in range(len(constraints)):
            if basis_vars[i] == -1:  # Ако няма базисна променлива в това ограничение
                artificial_var_index = num_vars + len(artificial_vars)
                for j in range(len(constraints)):
                    constraints[j][0].append(Rational(1) if j == i else Rational(0))  # Добавяне на изкуствена променлива
                artificial_vars.append(artificial_var_index)
                if objective == "min":
                    coefficients.append(Rational(M_VALUE))
                elif objective == "max":
                    coefficients.append(Rational(-M_VALUE)) # Изкуствените променливи, които се добавят, за да се създаде базис, влияят на целевата функция
    
        num_vars += len(artificial_vars)
    
        # Подготовка и извеждане на каноничния вид
        print("Каноничен вид на задачата:")
        print("Целева функция:")
        print_function(coefficients, num_vars, objective, free_coefficient)
    
        print("Ограничения:")
        for coeffs, sign, rhs in constraints:
            print_constraint(coeffs, rhs)

def canonicalize_for_table(coefficients, constraints, objective):
    num_vars = len(coefficients)
    surplus_vars = 0
    slack_vars = 0
    artificial_vars = 0

    # Добавяне на излишъци и изкуствени променливи
    for i, (coeffs, sign, rhs) in enumerate(constraints):
        if sign == "<=":
            coeffs.extend([Rational(1) if j == slack_vars else Rational(0) for j in range(len(constraints))])
            slack_vars += 1
        elif sign == ">=":
            coeffs.extend([Rational(-1) if j == surplus_vars else Rational(0) for j in range(len(constraints))])
            surplus_vars += 1
        constraints[i] = (coeffs, "=", rhs)

    # Проверка за базисни променливи
    basis_vars = find_basis(constraints)
    artificial_indices = []

    # Добавяне на изкуствени променливи при нужда
    for i in range(len(constraints)):
        if basis_vars[i] == -1:  # Ако няма базисна променлива в това ограничение
            artificial_var_index = num_vars + slack_vars + surplus_vars + artificial_vars
            for j in range(len(constraints)):
                constraints[j][0].append(Rational(1) if j == i else Rational(0))  # Добавяне на изкуствена променлива
            artificial_indices.append(artificial_var_index)
            artificial_vars += 1

    num_slack = slack_vars + surplus_vars
    total_vars = num_vars + num_slack + artificial_vars

    # Добавяне на ред с коефициентите от целевата функция
    objective_row = ['', '', ''] + coefficients + [Rational(0)] * (num_slack + artificial_vars)
    table = [objective_row]

    # Заглавия на колоните
    headers = ['CBx', 'Bx', 'b'] + [f"x{i + 1}" for i in range(total_vars)]
    table.append(headers)

    # Генериране на таблицата
    for i, (coeffs, sign, rhs) in enumerate(constraints):
        #base_var = f"x{total_vars - artificial_vars - i}" if surplus_vars > 0 else f"x{num_vars + i + 1 - len(constraints)}"
        #CBx = 0  # Коефициентът на базисната променлива в целевата функция
        basis_index = basis_vars[i] if basis_vars[i] != -1 else total_vars - artificial_vars - i
        base_var = f"x{basis_index + 1}"
        CBx = coefficients[basis_index] if basis_index < len(coefficients) else Rational(0)  # Взимаме коефициента от целевата функция
        bx_row = [CBx, base_var, Rational(rhs)] + coeffs
        table.append(bx_row)

    return table

def print_function(coeffs, num_vars, objective, free_coeff):
    terms = []
    for i, coeff in enumerate(coeffs[:num_vars]):
        if coeff == M_VALUE:
            terms.append(f"{str(coeff).replace(str(M_VALUE), 'M')}*x{i + 1}")
        elif coeff == -M_VALUE:
            terms.append(f"{str(coeff).replace(str(-M_VALUE), '-M')}*x{i + 1}")
        else:
            terms.append(f"{coeff}*x{i + 1}")
    
    # Присъединяване на свободния член само ако е различен от 0
    if free_coeff != 0:
        free_coeff_str = f" + {free_coeff}" if free_coeff > 0 else f" - {-free_coeff}"
    else:
        free_coeff_str = ""
    
    print(f"{'Максимизирай' if objective == 'max' else 'Минимизирай'}: {' + '.join(terms).replace('+ -', '- ')}", free_coeff_str)

def print_constraint(coeffs, rhs):
    terms = []
    for i, coeff in enumerate(coeffs):
        if coeff > 0:
            terms.append(f"{coeff}*x{i + 1}")
        elif coeff < 0:
            terms.append(f"- {abs(coeff)}*x{i + 1}")
        else:
            terms.append(f"0*x{i + 1}")

    # Заместване на първоначалния знак `+` с празен низ, ако е необходимо
    constraint_str = ' + '.join(terms).replace('+ -', '-')
    print(f"{constraint_str} = {rhs}")

def print_table(table, coefficients):
    approximation = simplify(sum(simplify(row[2]) * simplify(row[0].subs(M_VALUE, M)) for row in table[2:] if isinstance(row[0], (int, Rational, type(M)))))
    table[-1][2] = approximation
    for j in range(3, 3 + len(coefficients)):
        delta = simplify(sum(simplify(row[0].subs(M, M_VALUE)) * simplify(row[j]) for row in table[2:] if isinstance(row[0], (int, Rational, type(M)))) - coefficients[j - 3])
        if delta != 0:
            delta = simplify(sum(simplify(row[j] * simplify(row[0].subs(M_VALUE, M))) for row in table[2:] if isinstance(row[0], (int, Rational, type(M)))) - coefficients[j - 3].subs(M_VALUE, M))
        table[-1][j] = delta
    for row in table:
        print(" | ".join(f"{str(item).replace(str(M_VALUE), 'M'):>15}" for item in row))

def calculate_first_approximation(table, coefficients):
    """
    Пресмята първото приближение за целевата функция и делтите за всички небазисни променливи.
    """
    # Изчисляване на стойността на целевата функция Z
    numeric_coeffs = [c.subs(M, M_VALUE) if c.has(M) else c for c in coefficients]
    approximation = simplify(sum(simplify(row[0]).subs(M, M_VALUE) * simplify(row[2]) for row in table[2:] if isinstance(row[0], (int, Rational))))
    
    # Пресмятане на делтите
    deltas = []
    num_vars = len(coefficients)
    for j in range(3, 3 + num_vars):
        delta = simplify(sum(simplify(row[0]).subs(M, M_VALUE) * simplify(row[j]) for row in table[2:] if isinstance(row[0], (int, Rational))) - numeric_coeffs[j - 3])
        deltas.append(delta)

    return approximation, deltas

def add_approximation_and_deltas(table, coefficients):
    approximation, deltas = calculate_first_approximation(table, coefficients)
    approximation_row = [''] * 3 + deltas
    approximation_row[0] = "Z"
    approximation_row[1] = " = "
    approximation_row[2] = f" {approximation}"
    table.append(approximation_row)

def check_optimality(table, objective):
    # Взимаме делтите от последния ред на таблицата, изключвайки първите три елемента
    deltas = table[-1][3:]

    # За максимизация, проверяваме дали всички делти са >= 0
    if objective == "max":
        if all(delta.subs(M, M_VALUE) >= 0 for delta in deltas):
            print("\nРешението е оптимално.")
            return True
        else:
            print("\nРешението не е оптимално. Необходимо е допълнително итериране.")
            return False
    
    # За минимизация, проверяваме дали всички делти са <= 0
    elif objective == "min":
        if all(delta.subs(M, M_VALUE) <= 0 for delta in deltas):
            print("\nРешението е оптимално.")
            return True
        else:
            print("\nРешението не е оптимално. Необходимо е допълнително итериране.")
            return False

def find_pivot_element(table, objective):
    # Извличане на реда на делтите
    delta_row = table[-1][3:]  # Премахваме първите три колони ('CBx', 'Bx', 'b')
    for i in range(len(delta_row)):
        delta_row[i] = delta_row[i].subs(M, M_VALUE)
    
    # Определяне на най-големият нарушител на критерия за оптималност
    if objective == "max":
        # За максимизация, търсим най-малката (най-голяма отрицателна) делта
        pivot_col = min(enumerate(delta_row), key=lambda x: x[1])
    else:
        # За минимизация, търсим най-голямата (най-голяма положителна) делта
        pivot_col = max(enumerate(delta_row), key=lambda x: x[1])
    
    # Извличане на стойностите от стълба на избраната делта и стълба 'b'
    b_values = [row[2].subs(M, M_VALUE) for row in table[2:-1]]  # Вземаме всички стойности от стълба 'b', пропускайки заглавния ред и реда на делтите
    column_values = [row[pivot_col[0] + 3] for row in table[2:-1]]  # +3 за да компенсираме пропуснатите първите три колони
    
    if all(val <= 0 for val in column_values):
        print("Задачата няма решение, поради липса на положителен елемент в ключовия стълб")
        return None
    
    # Изчисляване на минималната дроб и определяне на ключовия елемент
    ratios = []
    for b, column_val in zip(b_values, column_values):
        if column_val > 0:  # Проверка за предотвратяване на деление на нула
            ratios.append(b / column_val)
        else:
            ratios.append(float('inf'))  # Неизползваемо висока стойност
    
    pivot_row_index = ratios.index(min(ratios))
    pivot_value = table[pivot_row_index + 2][pivot_col[0] + 3]  # +1 за да компенсираме пропуснатия заглавен ред
    
    return pivot_row_index + 2, pivot_col[0] + 3, pivot_value

# Забележка: Тази функция предполага, че 'table' вече съдържа канонизираната форма с добавените ред за делтите.
# Трябва да се адаптира спрямо точната структура на вашата таблица.

def pivot_table(table, pivot_row_index, pivot_col_index, coefficients):
    pivot_element = table[pivot_row_index][pivot_col_index]
    
    # Нормализиране на ключовия ред
    for i in range(2, len(table[pivot_row_index])):
        table[pivot_row_index][i] /= pivot_element
    
    # Обновяване на останалите редове
    for r in range(2, len(table) - 1):
        if r != pivot_row_index:
            factor = table[r][pivot_col_index]
            for c in range(2, len(table[r])):
                table[r][c] = table[r][c] - factor * table[pivot_row_index][c]
    
    # Заместване на променливата в Bx и обновяване на CBx
    bx_var = f"x{pivot_col_index - 2}"  # Адаптирайте според точната структура на таблицата
    table[pivot_row_index][1] = bx_var  # Заместваме в Bx
    table[pivot_row_index][0] = coefficients[pivot_col_index - 3]  # Обновяваме CBx със съответния коефициент от целевата функция

    return table

def print_optimal_solution(table, number_original_variables, objective, free_coeff):
    # Изчисляване на Z чрез обхождане на първи и трети стълб (индекси 0 и 2)
    z_value = Rational(0)
    
    for row in table[2:-1]:  # Пропускаме заглавния и последния ред
        z_value += simplify(row[0].subs(M_VALUE, M) * Rational(row[2]))  # Сумираме произведението на първи и трети стълб
    
    z_value += free_coeff
    
    contains_M_in_first_column = any(row[0].has(M_VALUE) for row in table[2:-1])
    if not z_value.has(M) and not contains_M_in_first_column:
        # Създаване на речник за оптималните стойности на променливите
        optimal_values = {f"x{var + 1}": Rational(0) for var in range(number_original_variables)}
        
        # Обхождане на редовете в таблицата за определяне на стойностите на променливите
        for row in table[1:-1]:  # Пропускаме първия (заглавния) и последния (делтите) ред
            bx_variable = row[1]  # Променливата в стълб Bx
            if bx_variable in optimal_values:
                optimal_values[bx_variable] = row[2]  # Стойността от стълб b
        
        # Извеждане на резултатите
        print(f"{objective.upper()} Z = {z_value}")
        for var, value in optimal_values.items():
            print(f"{var} = {value}", end=" ")
        print()
    else:
        # Създаване на речник за оптималните стойности на променливите
        optimal_values = {f"x{var - 2}": Rational(0) for var in range(3, len(table[1]))}
        
        # Обхождане на редовете в таблицата за определяне на стойностите на променливите
        for row in table[1:-1]:  # Пропускаме първия (заглавния) и последния (делтите) ред
            bx_variable = row[1]  # Променливата в стълб Bx
            if bx_variable in optimal_values:
                optimal_values[bx_variable] = row[2]  # Стойността от стълб b
        
        # Извеждане на резултатите
        print("\nРешение на М-задачата:")
        print(f"{objective.upper()} Z = {z_value}")
        for var, value in optimal_values.items():
            print(f"{var} = {value}", end=" ")
        print()
        print("Изходната задача няма решение, заради наличието на изкуствени променливи в базиса")
        
def main():
    num_vars = int(input("Въведете броя на променливите в задачата: "))
    coefficients = read_coefficients(num_vars, "целевата функция")
    
    free_coefficient = Rational(input("Въведете свободния коефициент, ако няма такъв, то въведете 0: "))
    
    objective = input("Въведете типа екстремум, който търсите (max за максимум или min за минимум): ")
    
    num_constraints = int(input("Въведете броя на ограниченията: "))
    constraints = read_constraints(num_vars, num_constraints)
    
    print("\nЦелева функция:")
    print(display_function(coefficients, objective, free_coefficient))
    
    print("\nОграничения:")
    print(display_constraints(constraints))
    
    print("\n")
    
    canonicalize(coefficients, constraints, objective, free_coefficient)
    canonical_form = canonicalize_for_table(coefficients, constraints, objective)
    print('\nПърва симплекс таблица\n')
    add_approximation_and_deltas(canonical_form, coefficients)
    print_table(canonical_form, coefficients)
    #check_optimality(canonical_form, objective)
    while check_optimality(canonical_form, objective) == False:
        if find_pivot_element(canonical_form, objective) is None:
            return
        pivot_row, pivot_col, pivot_value = find_pivot_element(canonical_form, objective)
        print(f"\nКлючовият елемент е: (ред: {canonical_form[pivot_row][1]}, стълб: {canonical_form[1][pivot_col]}) = {pivot_value}\n")
        print("Следваща симплекс таблица\n")
        canonical_form = pivot_table(canonical_form, pivot_row, pivot_col, coefficients)
        add_approximation_and_deltas(canonical_form, coefficients)
        del canonical_form[-2]
        print_table(canonical_form, coefficients)
    print_optimal_solution(canonical_form, num_vars, objective, free_coefficient)

if __name__ == "__main__":
    main()
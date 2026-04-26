from flask import Flask, render_template, request, jsonify

app = Flask(__name__)

# ==========================================
# Logic — Gaussian Elimination
# ==========================================

def forward_elimination(matrix):
    n = len(matrix)
    steps = []

    steps.append({
        "label": "Matrix เริ่มต้น [A|B]",
        "op": "",
        "matrix": [row[:] for row in matrix]
    })

    for pivot in range(n):
        if abs(matrix[pivot][pivot]) < 1e-10:
            for swap in range(pivot + 1, n):
                if abs(matrix[swap][pivot]) > 1e-10:
                    matrix[pivot], matrix[swap] = matrix[swap], matrix[pivot]
                    steps.append({
                        "label": f"สลับ R{pivot+1} กับ R{swap+1}",
                        "op": f"↕ R{pivot+1} ↔ R{swap+1}",
                        "matrix": [row[:] for row in matrix]
                    })
                    break

        for i in range(pivot + 1, n):
            if abs(matrix[pivot][pivot]) < 1e-10:
                continue
            factor = matrix[i][pivot] / matrix[pivot][pivot]
            if abs(factor) < 1e-10:
                continue
            for j in range(len(matrix[0])):
                matrix[i][j] -= factor * matrix[pivot][j]
            steps.append({
                "label": f"กำจัดคอลัมน์ {pivot+1}",
                "op": f"R{i+1} = R{i+1} - ({factor:.2f}) × R{pivot+1}",
                "matrix": [row[:] for row in matrix]
            })

    return matrix, steps


def check_solution(matrix):
    n = len(matrix)
    for i in range(n):
        all_zero_left = all(abs(matrix[i][j]) < 1e-10 for j in range(n))
        right_side = abs(matrix[i][-1]) < 1e-10
        if all_zero_left and not right_side:
            return "no_solution"
        if all_zero_left and right_side:
            return "infinite_solution"
    return "one_solution"


def back_substitution(matrix):
    n = len(matrix)
    solution = [0] * n
    for i in range(n - 1, -1, -1):
        solution[i] = matrix[i][-1]
        for j in range(i + 1, n):
            solution[i] -= matrix[i][j] * solution[j]
        solution[i] /= matrix[i][i]
    return solution


# ==========================================
# Logic — Inverse Matrix
# ==========================================

def get_determinant(A):
    n = len(A)
    if n == 1:
        return A[0][0]
    if n == 2:
        return A[0][0]*A[1][1] - A[0][1]*A[1][0]
    det = 0
    for c in range(n):
        sub = [row[:c] + row[c+1:] for row in A[1:]]
        det += ((-1) ** c) * A[0][c] * get_determinant(sub)
    return det


def get_cofactor_matrix(A):
    n = len(A)
    cofactor = []
    for r in range(n):
        row = []
        for c in range(n):
            sub = [A[i][:c] + A[i][c+1:] for i in range(n) if i != r]
            row.append(((-1) ** (r + c)) * get_determinant(sub))
        cofactor.append(row)
    return cofactor


def transpose(A):
    n = len(A)
    return [[A[j][i] for j in range(n)] for i in range(n)]


def mat_vec_multiply(A, b):
    return [sum(A[i][j] * b[j] for j in range(len(b))) for i in range(len(A))]


def solve_by_inverse(A, B):
    det = get_determinant(A)
    if abs(det) < 1e-10:
        return None, det, None, None, None

    cofactor = get_cofactor_matrix(A)
    adj = transpose(cofactor)
    inv = [[adj[i][j] / det for j in range(len(A))] for i in range(len(A))]
    solution = mat_vec_multiply(inv, B)

    return solution, det, cofactor, adj, inv


# ==========================================
# Routes
# ==========================================

@app.route('/')
def index():
    return render_template('index.html')


@app.route('/solve', methods=['POST'])
def solve():
    data = request.json
    raw = data['matrix']  # [[a,b,c,d], [a,b,c,d], [a,b,c,d]]
    n = len(raw)

    matrix = [row[:] for row in raw]
    A_orig = [row[:-1] for row in raw]
    B_orig = [row[-1] for row in raw]

    # Gaussian
    matrix, steps = forward_elimination(matrix)
    status = check_solution(matrix)

    result = {"status": status, "steps": steps}

    if status == "one_solution":
        sol_g = back_substitution(matrix)
        result["gaussian"] = sol_g

        # Inverse
        sol_i, det, cofactor, adj, inv = solve_by_inverse(A_orig, B_orig)
        if sol_i:
            result["inverse"] = {
                "solution": sol_i,
                "det": det,
                "cofactor": cofactor,
                "adjugate": adj,
                "inv": inv
            }

    return jsonify(result)


if __name__ == '__main__':
    app.run(debug=True)
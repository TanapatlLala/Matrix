from flask import Flask, render_template, request, jsonify  # นำเข้า Flask และฟังก์ชันที่ใช้รับ-ส่งข้อมูล

app = Flask(__name__)  # สร้างแอปพลิเคชัน Flask

# ==========================================
# Logic — Gaussian Elimination
# ==========================================

def forward_elimination(matrix):
    # ฟังก์ชันทำ Forward Elimination เพื่อแปลงเมทริกซ์เป็นรูป Row Echelon Form
    n = len(matrix)  # จำนวนแถว (= จำนวนตัวแปร)
    steps = []  # รายการเก็บขั้นตอนทั้งหมดเพื่อแสดงผล

    steps.append({
        "label": "Matrix เริ่มต้น [A|B]",  # ป้ายกำกับขั้นตอนแรก
        "op": "",                            # ยังไม่มีการดำเนินการ
        "matrix": [row[:] for row in matrix] # คัดลอกเมทริกซ์ปัจจุบันเก็บไว้
    })

    for pivot in range(n):  # วนซ้ำทุก pivot (แนวทแยงหลัก)

        if abs(matrix[pivot][pivot]) < 1e-10:  # ถ้า pivot เกือบเป็นศูนย์ (หารไม่ได้)
            for swap in range(pivot + 1, n):    # มองหาแถวด้านล่างที่ค่าไม่เป็นศูนย์
                if abs(matrix[swap][pivot]) > 1e-10:  # พบแถวที่ใช้แทนได้
                    matrix[pivot], matrix[swap] = matrix[swap], matrix[pivot]  # สลับแถว
                    steps.append({
                        "label": f"สลับ R{pivot+1} กับ R{swap+1}",  # บันทึกคำอธิบาย
                        "op": f"↕ R{pivot+1} ↔ R{swap+1}",          # บันทึกสัญลักษณ์การดำเนินการ
                        "matrix": [row[:] for row in matrix]          # บันทึกเมทริกซ์หลังสลับ
                    })
                    break  # สลับแล้วออกจากลูปหาแถว

        for i in range(pivot + 1, n):  # วนทุกแถวด้านล่าง pivot
            if abs(matrix[pivot][pivot]) < 1e-10:  # ถ้า pivot ยังเป็นศูนย์อยู่ ข้ามไป
                continue
            factor = matrix[i][pivot] / matrix[pivot][pivot]  # คำนวณตัวคูณ (multiplier)
            if abs(factor) < 1e-10:  # ถ้าตัวคูณเกือบศูนย์ ไม่ต้องทำอะไร
                continue
            for j in range(len(matrix[0])):  # วนทุกคอลัมน์ในแถวนั้น
                matrix[i][j] -= factor * matrix[pivot][j]  # ลบแถว pivot คูณด้วย factor ออก
            steps.append({
                "label": f"กำจัดคอลัมน์ {pivot+1}",                              # คำอธิบายขั้นตอน
                "op": f"R{i+1} = R{i+1} - ({factor:.2f}) × R{pivot+1}",           # สูตรที่ใช้
                "matrix": [row[:] for row in matrix]                               # เมทริกซ์หลังกำจัด
            })

    return matrix, steps  # คืนเมทริกซ์รูป Row Echelon และรายการขั้นตอน


def check_solution(matrix):
    # ตรวจสอบว่าระบบสมการมีคำตอบกี่แบบ
    n = len(matrix)  # จำนวนตัวแปร
    has_infinite = False  # ธงบอกว่ามีคำตอบอนันต์หรือไม่
    for i in range(n):  # ตรวจทุกแถว
        all_zero_left = all(abs(matrix[i][j]) < 1e-10 for j in range(n))  # ทุกค่าฝั่ง A เป็นศูนย์?
        right_side = abs(matrix[i][-1]) < 1e-10                            # ค่าฝั่ง B เป็นศูนย์?
        if all_zero_left and not right_side:  # แถวรูป [0 0 ... 0 | c≠0] → ขัดแย้ง
            return "no_solution"              # ไม่มีคำตอบ
        if all_zero_left and right_side:      # แถวรูป [0 0 ... 0 | 0] → ตัวแปรอิสระ
            has_infinite = True               # มีแนวโน้มคำตอบอนันต์
            
    if has_infinite:          # ถ้ามีแถวศูนย์ทั้งหมด
        return "infinite_solution"  # คำตอบอนันต์
        
    return "one_solution"  # มีคำตอบเดียว


def back_substitution(matrix):
    # แทนค่ากลับ (Back Substitution) เพื่อหาคำตอบจากเมทริกซ์รูป Row Echelon
    n = len(matrix)          # จำนวนตัวแปร
    solution = [0] * n       # ตัวแปรเก็บคำตอบ เริ่มต้นเป็นศูนย์ทั้งหมด
    back_steps = []          # รายการขั้นตอนการแทนค่ากลับ
    for i in range(n - 1, -1, -1):  # วนจากแถวล่างสุดขึ้นบน
        val = matrix[i][-1]          # ค่า RHS ของแถวนี้
        
        terms = []  # รายการเทอมที่ต้องนำไปลบ (จากตัวแปรที่รู้แล้ว)
        for j in range(i + 1, n):   # วนตัวแปรที่คำนวณไปแล้ว
            terms.append(f"({matrix[i][j]:.2f})*({solution[j]:.2f})")  # สร้างสตริงแสดงเทอม
        
        sub_str = " - ".join(terms) if terms else ""  # รวมเทอมเป็นสตริง
        if sub_str:  # ถ้ามีเทอมที่ต้องลบ
            step_str = f"x_{i+1} = ({val:.2f} - ({sub_str})) / {matrix[i][i]:.2f}"  # สร้างสูตรแสดง
        else:        # แถวล่างสุด ไม่มีเทอมอื่น
            step_str = f"x_{i+1} = {val:.2f} / {matrix[i][i]:.2f}"

        for j in range(i + 1, n):       # ลบค่าที่รู้แล้วออกจาก val จริงๆ
            val -= matrix[i][j] * solution[j]
        val /= matrix[i][i]              # หารด้วย coefficient ของตัวแปรตัวนี้
        solution[i] = val                # เก็บคำตอบ
        
        step_str += f" = {val:.4f}"      # เพิ่มผลลัพธ์ต่อท้ายสตริง
        back_steps.append(step_str)      # บันทึกขั้นตอน
        
    return solution, back_steps  # คืนคำตอบและรายการขั้นตอน


# ==========================================
# Logic — Inverse Matrix
# ==========================================

def get_determinant(A):
    # คำนวณ Determinant ของเมทริกซ์ A แบบ Recursive (Cofactor Expansion)
    n = len(A)
    if n == 1:  # กรณีฐาน: เมทริกซ์ 1×1
        return A[0][0]
    if n == 2:  # กรณีฐาน: เมทริกซ์ 2×2 → ad - bc
        return A[0][0]*A[1][1] - A[0][1]*A[1][0]
    det = 0  # เริ่ม det เป็นศูนย์
    for c in range(n):  # Cofactor Expansion ตามแถวแรก
        sub = [row[:c] + row[c+1:] for row in A[1:]]  # ตัดแถวแรกและคอลัมน์ c ออก (Minor)
        det += ((-1) ** c) * A[0][c] * get_determinant(sub)  # สะสม det ด้วยสูตร Cofactor
    return det  # คืนค่า determinant


def get_cofactor_matrix(A):
    # สร้าง Cofactor Matrix ของ A (ใช้ในการคำนวณ Adjugate และ Inverse)
    n = len(A)
    cofactor = []  # เก็บ cofactor ทั้งหมด
    for r in range(n):       # วนทุกแถว
        row = []
        for c in range(n):   # วนทุกคอลัมน์
            sub = [A[i][:c] + A[i][c+1:] for i in range(n) if i != r]  # Minor: ตัดแถว r และคอลัมน์ c
            row.append(((-1) ** (r + c)) * get_determinant(sub))        # คูณเครื่องหมาย (-1)^(r+c)
        cofactor.append(row)
    return cofactor  # คืน Cofactor Matrix


def transpose(A):
    # Transpose เมทริกซ์ (สลับแถวกับคอลัมน์)
    n = len(A)
    return [[A[j][i] for j in range(n)] for i in range(n)]  # สร้างเมทริกซ์ใหม่ที่สลับ index


def mat_vec_multiply(A, b):
    # คูณเมทริกซ์ A กับเวกเตอร์ b → คืนเวกเตอร์ผลลัพธ์
    return [sum(A[i][j] * b[j] for j in range(len(b))) for i in range(len(A))]


def build_det_step(A, det):
    """
    สร้างสตริงแสดงการคำนวณ det(A) แบบ "คูณลง - คูณขึ้น"
    - n=2 : สูตร ad - bc ตรงๆ
    - n=3 : Sarrus Rule (3 เส้นลง - 3 เส้นขึ้น)
    - n>3 : Cofactor Expansion พร้อมอธิบายทิศ
    """
    n = len(A)

    def fmt(v):
        # แสดงตัวเลขโดยไม่มีทศนิยมไม่จำเป็น
        return f"{v:.4g}"

    if n == 2:  # เมทริกซ์ 2×2 ใช้สูตร ad - bc โดยตรง
        a, b_ = A[0][0], A[0][1]  # แถวแรก
        c, d  = A[1][0], A[1][1]  # แถวที่สอง
        down = a * d   # ผลคูณเส้นลง (↘)
        up   = b_ * c  # ผลคูณเส้นขึ้น (↗)
        lines = [
            "┌ คูณลง (↘) ┐",
            f"  ({fmt(a)}) × ({fmt(d)}) = {fmt(down)}",
            "└ คูณขึ้น (↗) ┘",
            f"  ({fmt(b_)}) × ({fmt(c)}) = {fmt(up)}",
            "─────────────────",
            f"det(A) = {fmt(down)} − {fmt(up)} = {fmt(det)}"
        ]
        return "\n".join(lines)  # รวมบรรทัดเป็นสตริงเดียว

    elif n == 3:
        # Sarrus Rule: 3 เส้นลง (↘) หักด้วย 3 เส้นขึ้น (↗)
        a = A  # alias สั้นๆ ของเมทริกซ์
        d1 = a[0][0]*a[1][1]*a[2][2]  # เส้นทแยงหลัก
        d2 = a[0][1]*a[1][2]*a[2][0]  # เส้นที่สองลงขวา
        d3 = a[0][2]*a[1][0]*a[2][1]  # เส้นที่สามลงขวา
        u1 = a[0][2]*a[1][1]*a[2][0]  # เส้นทแยงรอง (ขึ้นซ้าย)
        u2 = a[0][0]*a[1][2]*a[2][1]  # เส้นที่สองขึ้นซ้าย
        u3 = a[0][1]*a[1][0]*a[2][2]  # เส้นที่สามขึ้นซ้าย
        sum_down = d1 + d2 + d3  # รวมผลคูณเส้นลง
        sum_up   = u1 + u2 + u3  # รวมผลคูณเส้นขึ้น
        lines = [
            "┌─── คูณลง (↘) ───────────────────────────────────────┐",
            f"  ({fmt(a[0][0])})×({fmt(a[1][1])})×({fmt(a[2][2])}) = {fmt(d1)}",
            f"  ({fmt(a[0][1])})×({fmt(a[1][2])})×({fmt(a[2][0])}) = {fmt(d2)}",
            f"  ({fmt(a[0][2])})×({fmt(a[1][0])})×({fmt(a[2][1])}) = {fmt(d3)}",
            f"  รวมคูณลง = {fmt(d1)} + {fmt(d2)} + {fmt(d3)} = {fmt(sum_down)}",
            "├─── คูณขึ้น (↗) ───────────────────────────────────────┤",
            f"  ({fmt(a[0][2])})×({fmt(a[1][1])})×({fmt(a[2][0])}) = {fmt(u1)}",
            f"  ({fmt(a[0][0])})×({fmt(a[1][2])})×({fmt(a[2][1])}) = {fmt(u2)}",
            f"  ({fmt(a[0][1])})×({fmt(a[1][0])})×({fmt(a[2][2])}) = {fmt(u3)}",
            f"  รวมคูณขึ้น = {fmt(u1)} + {fmt(u2)} + {fmt(u3)} = {fmt(sum_up)}",
            "├──────────────────────────────────────────────────────┤",
            f"  det(A) = {fmt(sum_down)} − {fmt(sum_up)} = {fmt(det)}"
        ]
        return "\n".join(lines)  # รวมเป็นสตริงเดียว

    else:
        # n>3: Cofactor Expansion แถวแรก พร้อมระบุทิศทาง
        cofactor = get_cofactor_matrix(A)  # คำนวณ Cofactor Matrix ก่อน
        lines = ["Cofactor Expansion (แถวแรก) — คูณลง/ขึ้น ตามเครื่องหมาย (+/−):"]
        terms = []  # เก็บเทอมแต่ละตัวเพื่อรวมท้ายสุด
        for c in range(n):
            sign = "+" if cofactor[0][c] >= 0 else "−"  # เครื่องหมายของ Cofactor
            sub = [A[i][:c] + A[i][c+1:] for i in range(1, n)]  # Minor matrix
            minor_val = get_determinant(sub)                      # det ของ Minor
            direction = "↘ (คูณลง)" if (c % 2 == 0) else "↗ (คูณขึ้น)"  # ทิศทางตาม index
            lines.append(
                f"  {direction}  a[1][{c+1}]×M{1}{c+1} = ({fmt(A[0][c])}) × ({fmt(minor_val)}) × ({'+1' if c%2==0 else '-1'}) = {fmt(A[0][c]*cofactor[0][c])}"
            )
            terms.append(fmt(A[0][c] * cofactor[0][c]))  # บันทึกค่าของเทอมนี้
        lines.append("─────────────────")
        lines.append(f"det(A) = {' + '.join(terms)} = {fmt(det)}")  # สรุปผล
        return "\n".join(lines)


def solve_by_inverse(A, B):
    # แก้สมการ Ax = B ด้วยวิธี Inverse Matrix: x = A⁻¹ × B
    n = len(A)
    det = get_determinant(A)             # คำนวณ determinant ของ A
    if abs(det) < 1e-10:                 # ถ้า det ≈ 0 แสดงว่าไม่มี inverse
        return None, det, None, None, None, None  # คืน None ทั้งหมด

    cofactor = get_cofactor_matrix(A)    # คำนวณ Cofactor Matrix
    det_step = build_det_step(A, det)    # สร้างข้อความแสดงขั้นตอนคำนวณ det

    adj = transpose(cofactor)            # Adjugate = Transpose ของ Cofactor
    inv = [[adj[i][j] / det for j in range(len(A))] for i in range(len(A))]  # A⁻¹ = adj / det
    solution = mat_vec_multiply(inv, B)  # x = A⁻¹ × B

    return solution, det, cofactor, adj, inv, det_step  # คืนคำตอบและทุกขั้นตอน


# ==========================================
# Routes
# ==========================================

@app.route('/')  # กำหนด URL หลัก (root) ของแอป
def index():
    return render_template('index.html')  # แสดงหน้า HTML หลัก


@app.route('/solve', methods=['POST'])  # รับ POST request ที่ /solve
def solve():
    data = request.json              # อ่าน JSON ที่ส่งมาจาก Frontend
    raw = data['matrix']             # ดึงเมทริกซ์ [A|B] จากข้อมูล
    n = len(raw)                     # จำนวนสมการ

    matrix = [row[:] for row in raw]         # คัดลอกเมทริกซ์สำหรับ Gaussian
    A_orig = [row[:-1] for row in raw]       # แยกส่วน A (ไม่รวม RHS) สำหรับ Inverse
    B_orig = [row[-1] for row in raw]        # แยกส่วน B (RHS เพียงอย่างเดียว)

    # ทำ Gaussian Elimination
    matrix, steps = forward_elimination(matrix)
    status = check_solution(matrix)          # ตรวจสอบประเภทคำตอบ

    result = {"status": status, "steps": steps}  # เตรียม dict ผลลัพธ์

    if status == "one_solution":  # ถ้ามีคำตอบเดียว
        sol_g, back_steps = back_substitution(matrix)  # คำนวณ Back Substitution
        result["gaussian"] = sol_g                      # เก็บคำตอบจาก Gaussian
        result["back_steps"] = back_steps               # เก็บขั้นตอน Back Substitution

        # คำนวณด้วย Inverse Matrix ด้วย
        sol_i, det, cofactor, adj, inv, det_step = solve_by_inverse(A_orig, B_orig)
        if sol_i:  # ถ้ามี Inverse (det ≠ 0)
            result["inverse"] = {
                "solution": sol_i,      # คำตอบจากวิธี Inverse
                "det": det,             # ค่า determinant
                "cofactor": cofactor,   # Cofactor Matrix
                "adjugate": adj,        # Adjugate Matrix
                "inv": inv,             # Inverse Matrix
                "det_step": det_step    # ขั้นตอนคำนวณ det แบบละเอียด
            }

    return jsonify(result)  # ส่งผลลัพธ์กลับเป็น JSON


if __name__ == '__main__':
    app.run(debug=True)  # รันแอปในโหมด debug (แสดง error โดยละเอียด)
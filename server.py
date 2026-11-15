"""SIMTELE-RISK · servidor local (solo librería estándar de Python) con SQLite.
Uso:  python3 server.py [puerto]   →  abrir http://127.0.0.1:8000
La base de datos queda en simtele.db (persiste al cerrar y reabrir)."""
import json, math, os, sqlite3, sys, time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import urlparse, parse_qs
BASE = os.path.dirname(os.path.abspath(__file__)); DB = os.path.join(BASE, "simtele.db")
def db():
    c = sqlite3.connect(DB); c.row_factory = sqlite3.Row; c.execute("PRAGMA foreign_keys=ON"); return c
def init():
    c = db(); c.executescript(open(os.path.join(BASE, "db", "schema.sql"), encoding="utf-8").read()); c.commit(); c.close()
def clean(o):
    if isinstance(o, float): return o if math.isfinite(o) else None
    if isinstance(o, dict): return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [clean(v) for v in o]
    return o
# ---------------- Campos, validación y valores predefinidos ----------------
# (clave, etiqueta, unidad, mín, máx, mín_estricto, opcional)
CAMPOS = {
 1: [("v0","Velocidad inicial v₀","m/s",0,100,True,False),("a","Desaceleración base a","m/s²",0,50,True,False),("b_arrastre","Término variable b (a(v)=a+b·v)","1/s",0,10,False,False),
     ("h","Paso h","s",0,5,True,False),("t_reac","Tiempo de reacción","s",0,10,False,False),("dist_disp","Distancia disponible (opcional)","m",0,10000,True,True)],
 2: [("D0","Deterioro inicial D₀","adim",0,1,False,False),("alpha","α (tasa lineal)","1/mes",0,1,False,False),("beta","β (tasa cuadrática)","1/mes²",0,0.1,False,False),
     ("Dumbral","Umbral de inspección D*","adim",0,1,True,False),("retraso","Retraso de mantenimiento","meses",0,120,False,False),
     ("tol","Tolerancia","",0,0.01,True,False),("maxit","Iteraciones máximas","",1,200000,False,False)],
 3: [("m","Masa m","kg",0,1e4,True,False),("c","Amortiguamiento c","N·s/m",0,1e4,False,False),("k","Rigidez k","N/m",0,1e6,True,False),
     ("x0","Posición inicial x₀","m",-100,100,False,False),("v0","Velocidad inicial v₀","m/s",-100,100,False,False),("F0","Amplitud de fuerza F₀","N",0,1e6,False,False),
     ("w","Frecuencia de fuerza ω","rad/s",0,100,False,False),("h","Paso h","s",0,0.1,True,False),("T","Tiempo total T","s",0,500,True,False)],
 4: [("a","Inicio del intervalo","",-100,100,False,False),("b","Fin del intervalo","",-100,100,False,False),("n","Número de nodos","",4,30,False,False),
     ("sigma","Ruido σ","",0,1,False,False),("seed","Semilla del ruido","",0,1e9,False,False),("ocultos","Puntos ocultos","",1,10,False,False)]}
def mat_norm(n): return [[4.0 if i==j else 1.0 if abs(i-j)==1 else 0.0 for j in range(n)] for i in range(n)]
def mat_est(n): return [[2.2 if i==j else -1.0 if abs(i-j)==1 else 0.0 for j in range(n)] for i in range(n)]
def mat_hil(n): return [[1/(i+j+1) for j in range(n)] for i in range(n)]
def mk(A, xs): return {"A": A, "b": [sum(a*x for a, x in zip(r, xs)) for r in A]}
PRESETS = {
 1: {"NORMAL": dict(v0=5,a=1,b_arrastre=0,h=0.05,t_reac=0,dist_disp=None), "ESTRES": dict(v0=6.5,a=0.6,b_arrastre=0.05,h=0.5,t_reac=1.5,dist_disp=12),
     "MAL_CONDICIONADO": dict(v0=5,a=1,b_arrastre=0,h=2.5,t_reac=0,dist_disp=None)},
 2: {"NORMAL": dict(D0=0.10,alpha=0.02,beta=0.001,Dumbral=1,retraso=0,tol=1e-10,maxit=20000, **mk(mat_norm(4),[1,2,3,4])),
     "ESTRES": dict(D0=0.10,alpha=0.022,beta=0.003,Dumbral=1,retraso=6,tol=1e-10,maxit=20000, **mk(mat_est(10),list(range(1,11)))),
     "MAL_CONDICIONADO": dict(D0=0.10,alpha=0.02,beta=0.001,Dumbral=1,retraso=0,tol=1e-10,maxit=200000, **mk(mat_hil(5),[1]*5))},
 3: {"NORMAL": dict(m=1,c=0.4,k=4,x0=0.1,v0=0,F0=0,w=2,h=0.01,T=20), "ESTRES": dict(m=1,c=0.2,k=4,x0=0.1,v0=0,F0=1,w=2,h=0.1,T=40),
     "MAL_CONDICIONADO": dict(m=1,c=4,k=4,x0=0.1,v0=0,F0=0,w=2,h=0.01,T=20)},
 4: {"NORMAL": dict(funcion="sin",a=0,b=6.283185,n=9,sigma=0,seed=1,ocultos=3), "ESTRES": dict(funcion="sin",a=0,b=6.283185,n=15,sigma=0.05,seed=1,ocultos=3),
     "MAL_CONDICIONADO": dict(funcion="exp",a=0,b=5,n=20,sigma=0.001,seed=1,ocultos=3)}}
def cholesky_ok(A):
    n = len(A); L = [[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1):
            s = A[i][j] - sum(L[i][k]*L[j][k] for k in range(j))
            if i == j:
                if s <= 0: return False
                L[i][i] = math.sqrt(s)
            else: L[i][j] = s / L[j][j]
    return True
def validar(caso, tipo, nombre, p):
    e = []
    if tipo not in ("NORMAL", "ESTRES", "MAL_CONDICIONADO"): e.append("Escenario inválido.")
    if not str(nombre or "").strip(): e.append("El nombre es obligatorio.")
    for k, lab, un, lo, hi, strict, opt in CAMPOS[caso]:
        v = p.get(k)
        if v in (None, ""):
            if opt: p[k] = None; continue
            e.append(f"{lab}: es obligatorio."); continue
        try: v = float(v)
        except (TypeError, ValueError): e.append(f"{lab}: debe ser un número."); continue
        if not math.isfinite(v): e.append(f"{lab}: valor no finito."); continue
        if (strict and v <= lo) or (not strict and v < lo) or v > hi:
            e.append(f"{lab}: debe estar en {'(' if strict else '['}{lo}, {hi}].")
        p[k] = v
    if caso == 4 and p.get("funcion") not in ("sin", "exp"): e.append("Función: elija sin o exp.")
    if e: return e
    if caso == 2:
        if p["Dumbral"] <= p["D0"]: e.append("El umbral D* debe ser mayor que D₀.")
        if p["alpha"] + p["beta"] <= 0: e.append("α y β no pueden ser ambos 0 (nunca se alcanzaría el umbral).")
        A, b = p.get("A"), p.get("b")
        try:
            A = [[float(x) for x in r] for r in A]; b = [float(x) for x in b]; n = len(b)
            if n < 4 or n > 30 or len(A) != n or any(len(r) != n for r in A): e.append("A debe ser cuadrada (4 ≤ n ≤ 30) y b del mismo tamaño.")
            elif any(abs(A[i][j]-A[j][i]) > 1e-9*max(1, abs(A[i][j])) for i in range(n) for j in range(n)): e.append("A debe ser simétrica (gradiente conjugado exige SPD).")
            elif not cholesky_ok(A): e.append("A no es definida positiva (SPD).")
            else: p["A"], p["b"] = A, b
        except Exception: e.append("A y b deben ser matrices/vectores numéricos.")
    if caso == 4 and p["b"] <= p["a"]: e.append("El fin del intervalo debe ser mayor que el inicio.")
    if caso == 4 and p["ocultos"] > p["n"] - 3: e.append("Demasiados puntos ocultos para ese número de nodos.")
    if caso == 3 and p["T"] / p["h"] > 2e5: e.append("T/h excede 200000 pasos.")
    if caso == 1 and p["h"] < 1e-5: e.append("h demasiado pequeño (mín. 1e-5).")
    return e
# ---------------- Métodos numéricos ----------------
def vd(a, b): return sum(x*y for x, y in zip(a, b))
def nr(a): return math.sqrt(vd(a, a))
def mv(A, x): return [vd(r, x) for r in A]
def gsolve(A, b):
    n = len(b); M = [r[:] + [b[i]] for i, r in enumerate(A)]
    for k in range(n):
        p = max(range(k, n), key=lambda i: abs(M[i][k])); M[k], M[p] = M[p], M[k]
        for i in range(k+1, n):
            f = M[i][k] / M[k][k]
            for j in range(k, n+1): M[i][j] -= f * M[k][j]
    x = [0.0]*n
    for i in range(n-1, -1, -1): x[i] = (M[i][n] - sum(M[i][j]*x[j] for j in range(i+1, n))) / M[i][i]
    return x
def eigsym(A):
    n = len(A); a = [r[:] for r in A]
    for _ in range(80):
        if sum(a[i][j]**2 for i in range(n) for j in range(i+1, n)) < 1e-28: break
        for p in range(n-1):
            for q in range(p+1, n):
                if abs(a[p][q]) < 1e-300: continue
                th = (a[q][q]-a[p][p]) / (2*a[p][q]); t = (1 if th >= 0 else -1) / (abs(th)+math.sqrt(th*th+1)); c = 1/math.sqrt(t*t+1); s = t*c
                for k in range(n): akp, akq = a[k][p], a[k][q]; a[k][p] = c*akp-s*akq; a[k][q] = s*akp+c*akq
                for k in range(n): apk, aqk = a[p][k], a[q][k]; a[p][k] = c*apk-s*aqk; a[q][k] = s*apk+c*aqk
    return sorted(a[i][i] for i in range(n))
def rres(A, x, b): return nr([bi - vd(r, x) for r, bi in zip(A, b)])
def iterativo(kind, A, b, tol, mx):
    n = len(b); x = [0.0]*n; nb = nr(b); R = []; X = []; ok = False
    if kind == "CG": r = b[:]; pv = r[:]; rs = vd(r, r)
    for k in range(1, mx+1):
        try:
            if kind == "JACOBI": x = [(b[i] - sum(A[i][j]*x[j] for j in range(n) if j != i)) / A[i][i] for i in range(n)]
            elif kind == "GAUSS_SEIDEL":
                for i in range(n): x[i] = (b[i] - sum(A[i][j]*x[j] for j in range(n) if j != i)) / A[i][i]
            else:
                Ap = mv(A, pv); al = rs / vd(pv, Ap); x = [x[i]+al*pv[i] for i in range(n)]; r = [r[i]-al*Ap[i] for i in range(n)]
                rn = vd(r, r); pv = [r[i] + (rn/rs)*pv[i] for i in range(n)]; rs = rn
            res = rres(A, x, b)
        except (OverflowError, ZeroDivisionError): break
        R.append(res); X.append(x[:] if k <= 300 else None)
        if not math.isfinite(res) or res > 1e100: break
        if res < tol * nb: ok = True; break
    return x, R, X, ok
def lu_solve(A, b):
    n = len(b); L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]; U = [r[:] for r in A]
    for k in range(n-1):
        for i in range(k+1, n):
            L[i][k] = U[i][k] / U[k][k]
            for j in range(k, n): U[i][j] -= L[i][k] * U[k][j]
    y = [0.0]*n
    for i in range(n): y[i] = b[i] - sum(L[i][j]*y[j] for j in range(i))
    x = [0.0]*n
    for i in range(n-1, -1, -1): x[i] = (y[i] - sum(U[i][j]*x[j] for j in range(i+1, n))) / U[i][i]
    return x
def errs(x, xr):
    ea = nr([a-b for a, b in zip(x, xr)]); n = nr(xr); return ea, (ea/n if n > 0 else None)
def relerr(e, ref): return e/abs(ref) if ref not in (0, None) and abs(ref) > 0 else None
def raices(p):
    D0, al, be, Du = p["D0"], p["alpha"], p["beta"], p["Dumbral"]; f = lambda t: D0 + al*t + be*t*t - Du; df = lambda t: al + 2*be*t
    tref = (Du-D0)/al if be < 1e-15 else (-al+math.sqrt(al*al+4*be*(Du-D0)))/(2*be)
    T = 1.0
    while f(T) < 0 and T < 1e7: T *= 2
    tol, mx = p["tol"], int(p["maxit"]); out = {}
    a_, b_ = 0.0, T; it = []
    for k in range(1, mx+1):
        c = (a_+b_)/2; it.append([c, abs(f(c))])
        if f(a_)*f(c) <= 0: b_ = c
        else: a_ = c
        if (b_-a_)/2 < tol: break
    out["BISECCION"] = (it, (b_-a_)/2 < tol, dict(a=0, b=T))
    x = T; it = []; okn = False
    for k in range(1, mx+1):
        d = df(x)
        if abs(d) < 1e-300: break
        x = x - f(x)/d; it.append([x, abs(f(x))])
        if abs(f(x)) < tol: okn = True; break
    out["NEWTON"] = (it, okn, dict(x0=T))
    x0, x1 = 0.0, T; it = []; oks = False
    for k in range(1, mx+1):
        f0, f1 = f(x0), f(x1)
        if f1 == f0: break
        x2 = x1 - f1*(x1-x0)/(f1-f0); x0, x1 = x1, x2; it.append([x2, abs(f(x2))])
        if abs(f(x2)) < tol: oks = True; break
    out["SECANTE"] = (it, oks, dict(x0=0, x1=T))
    return tref, out
def rk(order, f, y, t, h):
    add = lambda y, k, s: [yi + s*ki for yi, ki in zip(y, k)]
    if order == 2: k2 = f(t+h/2, add(y, f(t, y), h/2)); return add(y, k2, h)
    k1 = f(t, y); k2 = f(t+h/2, add(y, k1, h/2)); k3 = f(t+h/2, add(y, k2, h/2)); k4 = f(t+h, add(y, k3, h))
    return [y[i] + h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]) for i in range(len(y))]
def det_ref(p):
    v0, a, b = p["v0"], p["a"], p["b_arrastre"]
    if b < 1e-12: td = v0/a; xf = lambda t: v0*t - a*t*t/2
    else: td = math.log(1 + b*v0/a)/b; xf = lambda t: (v0 + a/b)*(1-math.exp(-b*t))/b - a*t/b
    return td, xf
def rk_run(order, p, h):
    a, b = p["a"], p["b_arrastre"]; td, xf = det_ref(p); f = lambda t, y: [y[1], -(a + b*y[1])]
    t = 0.0; y = [0.0, p["v0"]]; rows = [[0, 0.0, y[0], y[1], 0.0]]; k = 0
    while t < td - 1e-12:
        hh = min(h, td - t); y = rk(order, f, y, t, hh); t += hh; k += 1
        rows.append([k, t, y[0], y[1], abs(y[0]-xf(t))])
    return y, rows, xf, td
# ---------------- Ejecución y persistencia ----------------
def ins_ej(c, sid, metodo, tol, paso, mx, conv, ms):
    return c.execute("INSERT INTO ejecucion_metodo(simulacion_id,metodo,tolerancia,paso,max_iter,convergio,tiempo_ms) VALUES(?,?,?,?,?,?,?)", (sid, metodo, tol, paso, mx, int(conv), ms)).lastrowid
def ins_res(c, eid, ref, val, met, interp=None):
    c.execute("INSERT INTO resultado_resumen(ejecucion_id,referencia,valor_json,metrica_error_json,interpretacion) VALUES(?,?,?,?,?)", (eid, ref, json.dumps(clean(val)), json.dumps(clean(met)), interp))
def ins_it(c, eid, rows):
    c.executemany("INSERT INTO iteracion_resultado(ejecucion_id,numero,aproximacion_json,residual,error_abs,error_rel) VALUES(?,?,?,?,?,?)", [(eid, n, json.dumps(clean(a)), clean(r), clean(ea), clean(er)) for n, a, r, ea, er in rows])
def fmt(x): return "n/d" if x is None else f"{x:.3e}"
def ejecutar(c, sid, caso, p):
    concl = {"numerica": "", "problema": ""}; vals = {}
    if caso == 1:
        td, xf = det_ref(p); dref = xf(td); dtot = dref + p["v0"]*p["t_reac"]; ttot = td + p["t_reac"]; res = {}
        for order, nm in ((2, "RK2_PUNTO_MEDIO"), (4, "RK4")):
            t0 = time.time(); y, rows, _, _ = rk_run(order, p, p["h"]); ms = (time.time()-t0)*1000
            step = max(1, len(rows)//300); sel = [r for i, r in enumerate(rows) if i % step == 0 or i == len(rows)-1]
            eid = ins_ej(c, sid, nm, None, p["h"], len(rows)-1, True, ms)
            ins_it(c, eid, [(int(r[0]), [r[1], r[2], r[3]], None, r[4], relerr(r[4], xf(r[1]))) for r in sel])
            ex = abs(y[0]-dref); ev = abs(y[1]); tab = []
            for h in (p["h"], p["h"]/2, p["h"]/4, p["h"]/8):
                yy, _, _, _ = rk_run(order, p, h); tab.append([h, abs(yy[0]-dref)])
            ordp = [math.log2(tab[i][1]/tab[i+1][1]) if tab[i][1] > 1e-13 and tab[i+1][1] > 1e-13 else None for i in range(3)]
            ins_res(c, eid, f"analítica: t_det={td:.10g}, d_frenado={dref:.10g}", {"x_final": y[0], "v_final": y[1], "refinamiento": tab, "orden_estimado": ordp},
                    {"error_abs_x": ex, "error_rel_x": relerr(ex, dref), "error_abs_v": ev}, None)
            vals[nm] = ex; res[nm] = (ex, ordp)
        eid = ins_ej(c, sid, "REFERENCIA_ANALITICA", None, None, None, True, 0)
        marg = None if p.get("dist_disp") is None else p["dist_disp"] - dtot
        ins_res(c, eid, "v(t), x(t) en forma cerrada", {"t_det": td, "d_frenado": dref, "d_total_con_reaccion": dtot, "t_total": ttot, "margen": marg}, {}, None)
        concl["numerica"] = (f"Error en la distancia final: RK2 = {fmt(res['RK2_PUNTO_MEDIO'][0])}, RK4 = {fmt(res['RK4'][0])} (paso h = {p['h']}). "
            + ("Con desaceleración constante la solución es un polinomio de grado 2 y ambos métodos son prácticamente exactos; use b > 0 para ver la diferencia de orden. " if p["b_arrastre"] < 1e-12
               else f"Orden estimado al refinar h: RK2 ≈ {[round(o,2) for o in res['RK2_PUNTO_MEDIO'][1] if o]}, RK4 ≈ {[round(o,2) for o in res['RK4'][1] if o]}. ")
            + "Las integrales por trapecio, Simpson y Romberg y las raíces se añaden en la Semana 2.")
        concl["problema"] = (f"Con v₀ = {p['v0']} m/s el tiempo de detención es {td:.3f} s y la distancia de frenado {dref:.3f} m; con {p['t_reac']} s de reacción la distancia total es {dtot:.3f} m en {ttot:.3f} s. "
            + ("No se indicó distancia disponible. " if marg is None else (f"Quedan {marg:.2f} m de margen sobre la distancia disponible. " if marg >= 0 else f"Se excede la distancia disponible por {-marg:.2f} m. "))
            + "Modelo hipotético con datos sintéticos: no permite inferir el comportamiento de un sistema real.")
    elif caso == 2:
        tref, rt = raices(p); A, b = p["A"], p["b"]; f = lambda t: p["D0"] + p["alpha"]*t + p["beta"]*t*t
        txt = []
        for nm, (it, ok, extra) in rt.items():
            eid = ins_ej(c, sid, nm, p["tol"], None, int(p["maxit"]), ok, 0)
            ins_it(c, eid, [(i+1, [r[0]], r[1], abs(r[0]-tref), relerr(abs(r[0]-tref), tref)) for i, r in enumerate(it[:300])])
            last = it[-1][0] if it else None; ea = abs(last-tref) if last is not None else None
            ins_res(c, eid, f"analítica t*={tref:.12g} meses", {"t_estimado": last, "iteraciones": len(it), **extra}, {"error_abs": ea, "error_rel": relerr(ea, tref) if ea is not None else None}, None)
            vals[nm] = last; txt.append(f"{nm.capitalize()}: {len(it)} it., error {fmt(ea)}")
        lam = eigsym(A); kap = lam[-1]/lam[0]; xr = gsolve(A, b); lin = {}
        for kind in ("JACOBI", "GAUSS_SEIDEL", "CG", "LU"):
            t0 = time.time()
            if kind == "LU": x = lu_solve(A, b); R = [rres(A, x, b)]; X = [x]; ok = True
            else: x, R, X, ok = iterativo(kind, A, b, p["tol"], int(p["maxit"]))
            ms = (time.time()-t0)*1000; ea, er = errs(x, xr) if all(math.isfinite(v) for v in x) else (None, None)
            eid = ins_ej(c, sid, kind, p["tol"], None, int(p["maxit"]), ok, ms)
            rows = [(k+1, X[k], R[k], errs(X[k], xr)[0], errs(X[k], xr)[1]) for k in range(len(R)) if X[k] is not None and all(math.isfinite(v) for v in X[k])]
            if not rows: rows = [(0, [None], None, None, None)]
            ins_it(c, eid, rows)
            ins_res(c, eid, "solución directa (eliminación con pivoteo)", {"x": x if ok or kind == "LU" else None, "iteraciones": len(R), "kappa2": kap}, {"error_abs": ea, "error_rel": er, "residual": R[-1] if R else None}, None)
            lin[kind] = (ok, len(R), er, R[-1] if R else None); vals[kind] = er
        eid = ins_ej(c, sid, "REFERENCIA_ANALITICA", None, None, None, True, 0)
        dr = f(tref + p["retraso"]); sens = -tref/(p["alpha"] + 2*p["beta"]*tref) if (p["alpha"] + 2*p["beta"]*tref) else None
        ins_res(c, eid, "raíz positiva de la cuadrática", {"t_umbral": tref, "D_en_t_mas_retraso": dr, "dt_dalpha": sens, "kappa2": kap}, {}, None)
        nc = [k for k, v in lin.items() if not v[0]]
        concl["numerica"] = ("Raíces — " + "; ".join(txt) + ". Sistema lineal (n=%d, κ₂≈%.3g): " % (len(b), kap)
            + ", ".join(f"{k}: {'sí' if v[0] else 'NO'} converge, {v[1]} it." for k, v in lin.items()) + ". "
            + (f"No convergen: {', '.join(nc)} (fallo de convergencia, p. ej. ρ(B_J) ≥ 1, distinto del mal condicionamiento). " if nc else "")
            + ("κ elevado: el residuo pequeño no garantiza error pequeño; pérdida ≈ %.1f cifras." % math.log10(kap) if kap > 1e4 else "Sistema bien condicionado."))
        concl["problema"] = (f"El indicador alcanza el umbral D* = {p['Dumbral']} en t* ≈ {tref:.2f} meses. " + (f"Si el mantenimiento se retrasa {p['retraso']:g} meses, D(t*+retraso) = {dr:.3f} (excede el umbral). " if p["retraso"] > 0 else "")
            + (f"Sensibilidad: ∂t*/∂α ≈ {sens:.2f} meses por unidad de α. " if sens else "") + "Los coeficientes calibrados con un sistema mal condicionado propagan error a t*. Indicador hipotético: no es una predicción de fallas reales.")
    elif caso == 3:
        m, cc, k, x0, v0, F0, w, h, T = (p[q] for q in ("m", "c", "k", "x0", "v0", "F0", "w", "h", "T")); w0 = math.sqrt(k/m); z = cc/(2*math.sqrt(k*m)); wd = w0*math.sqrt(max(0, 1-z*z))
        reg = "subamortiguado" if z < 1-1e-9 else ("crítico" if abs(z-1) <= 1e-9 else "sobreamortiguado"); amp = F0/math.sqrt((k-m*w*w)**2 + (cc*w)**2) if F0 > 0 else None
        eid = ins_ej(c, sid, "REFERENCIA_ANALITICA", None, None, None, True, 0)
        ins_res(c, eid, "oscilador lineal de 1 g.d.l.", {"omega0": w0, "zeta": z, "omega_d": wd if z < 1 else None, "periodo_d": 2*math.pi/wd if 0 < wd else None, "regimen": reg, "amplitud_estacionaria": amp}, {}, None)
        concl["numerica"] = "Los métodos RK2/RK4 y la integración de x(t)² (trapecio, Simpson, Romberg) se ejecutan en la Semana 3. Por ahora se calculan los parámetros analíticos de referencia."
        concl["problema"] = (f"ω₀ = {w0:.4f} rad/s, ζ = {z:.4f} → régimen {reg}. " + (f"Amplitud estacionaria bajo F₀ sinusoidal: {amp:.4g} m. " if amp else "Sin fuerza externa (vibración libre). ")
            + ("ζ ≈ 1: raíces características casi iguales; sensibilidad alta a (c,k). " if abs(z-1) < 0.05 else "") + "Modelo lineal hipotético: no incluye no linealidades ni efectos reales del viento.")
    else: raise ValueError("La ejecución numérica del módulo Sensores (interpolación) se implementa en la Semana 2.")
    eid = ins_ej(c, sid, "CONCLUSION", None, None, None, True, 0); ins_res(c, eid, None, {}, {}, json.dumps(concl, ensure_ascii=False))
    return vals
def prev_vals(c, sid):
    out = {}
    for r in c.execute("SELECT e.metodo, r.valor_json, r.metrica_error_json FROM ejecucion_metodo e JOIN resultado_resumen r ON r.ejecucion_id=e.id WHERE e.simulacion_id=?", (sid,)):
        out[r["metodo"]] = (json.loads(r["valor_json"] or "{}"), json.loads(r["metrica_error_json"] or "{}"))
    return out
def scal(d):
    s = {}
    for k, v in d.items():
        if isinstance(v, (int, float)) and not isinstance(v, bool): s[k] = float(v)
    return s
def comparar(old, new):
    out = {}
    for m in new:
        if m in old and m != "CONCLUSION":
            a, b = scal(old[m][0]), scal(new[m][0]); keys = [k for k in a if k in b]
            out[m] = max((abs(a[k]-b[k]) for k in keys), default=0.0)
    return out
# ---------------- Acceso a datos ----------------
def problema_row(c, sid):
    r = c.execute("SELECT s.id,s.nombre,s.notas,s.fecha,s.estado,s.equipo,e.id AS esc_id,e.tipo,e.parametros_json,e.fuente_datos,e.caso_id,k.nombre AS caso FROM simulacion s JOIN escenario e ON e.id=s.escenario_id JOIN caso_estudio k ON k.id=e.caso_id WHERE s.id=?", (sid,)).fetchone()
    if not r: return None
    d = dict(r); d["parametros"] = json.loads(d.pop("parametros_json")); return d
def detalle(c, sid):
    d = problema_row(c, sid)
    if not d: return None
    ex = []; concl = None
    for e in c.execute("SELECT * FROM ejecucion_metodo WHERE simulacion_id=? ORDER BY id", (sid,)):
        e = dict(e); rs = c.execute("SELECT * FROM resultado_resumen WHERE ejecucion_id=?", (e["id"],)).fetchone()
        if e["metodo"] == "CONCLUSION": concl = json.loads(rs["interpretacion"]); continue
        e["resumen"] = {"referencia": rs["referencia"], "valor": json.loads(rs["valor_json"] or "{}"), "error": json.loads(rs["metrica_error_json"] or "{}")} if rs else None
        its = [dict(i) for i in c.execute("SELECT numero,aproximacion_json,residual,error_abs,error_rel FROM iteracion_resultado WHERE ejecucion_id=? ORDER BY numero", (e["id"],))]
        for i in its: i["aproximacion"] = json.loads(i.pop("aproximacion_json"))
        e["iteraciones"] = its; ex.append(e)
    d["ejecuciones"] = ex; d["conclusiones"] = concl; return d
class H(BaseHTTPRequestHandler):
    def log_message(self, *a): pass
    def send(self, code, obj, ctype="application/json"):
        body = (json.dumps(clean(obj), ensure_ascii=False) if ctype == "application/json" else obj).encode("utf-8")
        self.send_response(code); self.send_header("Content-Type", ctype + "; charset=utf-8"); self.send_header("Content-Length", str(len(body))); self.end_headers(); self.wfile.write(body)
    def body(self):
        n = int(self.headers.get("Content-Length") or 0); return json.loads(self.rfile.read(n) or b"{}")
    def handle_any(self, method):
        u = urlparse(self.path); q = parse_qs(u.query); parts = [x for x in u.path.split("/") if x]; c = db()
        try:
            if method == "GET" and not parts:
                return self.send(200, open(os.path.join(BASE, "index.html"), encoding="utf-8").read(), "text/html")
            if parts[:2] == ["api", "casos"] and method == "GET":
                rows = [dict(r) for r in c.execute("SELECT k.id,k.nombre,(SELECT COUNT(*) FROM escenario e JOIN simulacion s ON s.escenario_id=e.id WHERE e.caso_id=k.id) AS total FROM caso_estudio k ORDER BY k.id")]
                for r in rows: r["campos"] = [dict(k=a, label=b, unit=u_, min=lo, max=hi, strict=s, opt=o) for a, b, u_, lo, hi, s, o in CAMPOS[r["id"]]]; r["presets"] = PRESETS[r["id"]]
                return self.send(200, rows)
            if parts[:2] == ["api", "problemas"]:
                if len(parts) == 2 and method == "GET":
                    sql = "SELECT s.id,s.nombre,s.fecha,s.estado,e.tipo,e.caso_id,k.nombre AS caso FROM simulacion s JOIN escenario e ON e.id=s.escenario_id JOIN caso_estudio k ON k.id=e.caso_id WHERE 1=1"; a = []
                    if q.get("caso"): sql += " AND e.caso_id=?"; a.append(int(q["caso"][0]))
                    if q.get("q"): sql += " AND (s.nombre LIKE ? OR s.notas LIKE ? OR e.tipo LIKE ? OR k.nombre LIKE ?)"; a += ["%" + q["q"][0] + "%"] * 4
                    return self.send(200, [dict(r) for r in c.execute(sql + " ORDER BY s.id DESC", a)])
                if len(parts) == 2 and method == "POST":
                    d = self.body(); caso = int(d.get("caso_id", 0)); p = d.get("parametros") or {}
                    if caso not in CAMPOS: return self.send(400, {"errores": ["Caso inválido."]})
                    er = validar(caso, d.get("tipo"), d.get("nombre"), p)
                    if er: return self.send(400, {"errores": er})
                    eid = c.execute("INSERT INTO escenario(caso_id,tipo,parametros_json,fuente_datos) VALUES(?,?,?,?)", (caso, d["tipo"], json.dumps(p), d.get("fuente_datos") or "supuesto del equipo")).lastrowid
                    sid = c.execute("INSERT INTO simulacion(escenario_id,nombre,notas,equipo) VALUES(?,?,?,?)", (eid, d["nombre"].strip(), d.get("notas"), d.get("equipo"))).lastrowid
                    c.commit(); return self.send(201, problema_row(c, sid))
                sid = int(parts[2]); cur = problema_row(c, sid)
                if not cur: return self.send(404, {"errores": ["No existe el problema."]})
                if len(parts) == 3 and method == "GET": return self.send(200, detalle(c, sid))
                if len(parts) == 3 and method == "PUT":
                    d = self.body(); p = d.get("parametros") or {}; er = validar(cur["caso_id"], d.get("tipo"), d.get("nombre"), p)
                    if er: return self.send(400, {"errores": er})
                    c.execute("UPDATE escenario SET tipo=?,parametros_json=?,fuente_datos=? WHERE id=?", (d["tipo"], json.dumps(p), d.get("fuente_datos") or cur["fuente_datos"], cur["esc_id"]))
                    c.execute("UPDATE simulacion SET nombre=?,notas=?,estado='PENDIENTE' WHERE id=?", (d["nombre"].strip(), d.get("notas"), sid))
                    c.execute("DELETE FROM ejecucion_metodo WHERE simulacion_id=?", (sid,)); c.commit(); return self.send(200, problema_row(c, sid))
                if len(parts) == 3 and method == "DELETE":
                    c.execute("DELETE FROM escenario WHERE id=?", (cur["esc_id"],)); c.commit(); return self.send(200, {"eliminado": sid})
                if len(parts) == 4 and parts[3] == "ejecutar" and method == "POST":
                    old = prev_vals(c, sid)
                    c.execute("DELETE FROM ejecucion_metodo WHERE simulacion_id=?", (sid,))
                    try: ejecutar(c, sid, cur["caso_id"], cur["parametros"])
                    except ValueError as ex: c.rollback(); return self.send(501, {"errores": [str(ex)]})
                    c.execute("UPDATE simulacion SET estado='COMPLETA' WHERE id=?", (sid,)); c.commit()
                    out = detalle(c, sid); out["comparacion"] = comparar(old, prev_vals(c, sid)); return self.send(200, out)
            return self.send(404, {"errores": ["Ruta no encontrada."]})
        except Exception as ex:
            c.rollback(); return self.send(500, {"errores": [f"Error interno: {ex}"]})
        finally: c.close()
    def do_GET(self): self.handle_any("GET")
    def do_POST(self): self.handle_any("POST")
    def do_PUT(self): self.handle_any("PUT")
    def do_DELETE(self): self.handle_any("DELETE")
if __name__ == "__main__":
    init(); port = int(sys.argv[1]) if len(sys.argv) > 1 else 8000
    print(f"SIMTELE-RISK en http://127.0.0.1:{port}  (base de datos: {DB})"); ThreadingHTTPServer(("127.0.0.1", port), H).serve_forever()

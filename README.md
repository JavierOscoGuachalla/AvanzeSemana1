# SIMTELE-RISK

Desafío Integrador Final de **Métodos Numéricos**: aplicación web con cuatro casos de estudio inspirados en un sistema de transporte por cable (La Paz). Los escenarios son **hipotéticos** y no describen fallas verificadas de ninguna empresa.

## Módulos
1. **Detención de emergencia** — RK2/RK4, integración, raíces.
2. **Deterioro y mantenimiento** — raíces y sistemas lineales (Jacobi, Gauss-Seidel, LU, CG).
3. **Oscilaciones por viento** — RK2/RK4, integración.
4. **Reconstrucción de sensores** — Lagrange, Newton, splines cúbicos.

Cada módulo trabaja con los escenarios **NORMAL, ESTRÉS y MAL CONDICIONADO** y aplica teoría de errores (absoluto, relativo, porcentual, residual, sensibilidad).

## Ejecutar
Requiere solo Python 3 (sin dependencias externas).
```bash
python3 server.py            # crea simtele.db a partir de db/schema.sql y levanta el servidor
# abrir http://127.0.0.1:8000
```
Los datos se guardan en SQLite (`simtele.db`) y persisten al cerrar y reabrir el servidor.

## CRUD de problemas
| Operación | Endpoint |
|---|---|
| Crear | `POST /api/problemas` |
| Listar / buscar | `GET /api/problemas?caso=&q=` |
| Ver detalle | `GET /api/problemas/{id}` |
| Modificar | `PUT /api/problemas/{id}` |
| Eliminar | `DELETE /api/problemas/{id}` |
| Ejecutar / relanzar | `POST /api/problemas/{id}/ejecutar` |

## Estado (Semana 1)
- [x] Prototipo web con 4 módulos y CRUD con validaciones
- [x] Teoría de errores + sistemas lineales + raíces (bisección, Newton, secante) en el Módulo 2
- [x] RK2/RK4 y referencia analítica en el Módulo 1; parámetros analíticos en el Módulo 3
- [x] DER y SQL (`db/`), documento de plan y avance (`docs/SIMTELE-RISK_Plan_Semana1.docx`)
- [ ] Semana 2: interpolación (Módulo 4), integración numérica, raíces de v(t)
- [ ] Semana 3: RK en oscilaciones, integración final, prueba ejecutar → guardar → recuperar → repetir

## Estructura
```
server.py              servidor + API REST + métodos numéricos + SQLite
index.html             interfaz (inicio, módulos, formularios, detalle, historial)
db/schema.sql, DER.md  modelo de datos
python/                lineales.py, generar_avance.py (PDF), generar_word.js (Word)
docs/                  documento Word, PDF, figuras
```

## Equipo
Integrantes: _completar_

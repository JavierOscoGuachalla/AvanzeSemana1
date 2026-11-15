-- SIMTELE-RISK · Modelo relacional (SQLite). Se aplica automáticamente al iniciar server.py
CREATE TABLE IF NOT EXISTS caso_estudio (
  id             INTEGER PRIMARY KEY,
  nombre         TEXT NOT NULL UNIQUE,
  modelo_version TEXT NOT NULL DEFAULT '1.0'
);
CREATE TABLE IF NOT EXISTS escenario (
  id              INTEGER PRIMARY KEY,
  caso_id         INTEGER NOT NULL REFERENCES caso_estudio(id) ON DELETE CASCADE,
  tipo            TEXT NOT NULL CHECK (tipo IN ('NORMAL','ESTRES','MAL_CONDICIONADO')),
  parametros_json TEXT NOT NULL,
  fuente_datos    TEXT NOT NULL DEFAULT 'supuesto del equipo'
);
CREATE TABLE IF NOT EXISTS simulacion (
  id           INTEGER PRIMARY KEY,
  escenario_id INTEGER NOT NULL REFERENCES escenario(id) ON DELETE CASCADE,
  nombre       TEXT NOT NULL,
  notas        TEXT,
  fecha        TEXT NOT NULL DEFAULT CURRENT_TIMESTAMP,
  equipo       TEXT,
  estado       TEXT NOT NULL DEFAULT 'PENDIENTE' CHECK (estado IN ('PENDIENTE','COMPLETA','ERROR'))
);
CREATE TABLE IF NOT EXISTS ejecucion_metodo (
  id            INTEGER PRIMARY KEY,
  simulacion_id INTEGER NOT NULL REFERENCES simulacion(id) ON DELETE CASCADE,
  metodo        TEXT NOT NULL,
  tolerancia    REAL,
  paso          REAL,
  max_iter      INTEGER,
  convergio     INTEGER NOT NULL CHECK (convergio IN (0,1)),
  tiempo_ms     REAL
);
CREATE TABLE IF NOT EXISTS iteracion_resultado (
  id                INTEGER PRIMARY KEY,
  ejecucion_id      INTEGER NOT NULL REFERENCES ejecucion_metodo(id) ON DELETE CASCADE,
  numero            INTEGER NOT NULL,
  aproximacion_json TEXT NOT NULL,
  residual          REAL,
  error_abs         REAL,
  error_rel         REAL
);
CREATE TABLE IF NOT EXISTS resultado_resumen (
  id                 INTEGER PRIMARY KEY,
  ejecucion_id       INTEGER NOT NULL REFERENCES ejecucion_metodo(id) ON DELETE CASCADE,
  referencia         TEXT,
  valor_json         TEXT,
  metrica_error_json TEXT,
  interpretacion     TEXT
);
CREATE INDEX IF NOT EXISTS idx_esc_caso ON escenario(caso_id);
CREATE INDEX IF NOT EXISTS idx_sim_esc  ON simulacion(escenario_id);
CREATE INDEX IF NOT EXISTS idx_ej_sim   ON ejecucion_metodo(simulacion_id);
CREATE INDEX IF NOT EXISTS idx_it_ej    ON iteracion_resultado(ejecucion_id);
INSERT OR IGNORE INTO caso_estudio(id,nombre) VALUES (1,'Detención'),(2,'Deterioro'),(3,'Oscilaciones'),(4,'Sensores');

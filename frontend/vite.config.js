import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// launch.py (run.bat / run.sh) picks free ports from ../port.json and passes them in through these
// variables. Running `npm run dev` directly falls back to 5173 -> backend 8000.
const backendPort = Number(process.env.GENOMESCOPE_BACKEND_PORT) || 8000
const frontendPort = Number(process.env.GENOMESCOPE_FRONTEND_PORT) || 5173

export default defineConfig({
  plugins: [react()],
  server: {
    port: frontendPort,
    strictPort: true, // fail loudly instead of silently moving to a port the backend's CORS list doesn't know
    proxy: {
      '/api': {
        target: `http://127.0.0.1:${backendPort}`,
        changeOrigin: true,
      },
    },
  },
  preview: {
    port: frontendPort,
    strictPort: true,
  },
})

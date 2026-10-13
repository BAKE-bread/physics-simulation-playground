import { createServer } from 'node:http';
import { readFile, stat } from 'node:fs/promises';
import { extname, join, normalize, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = resolve(fileURLToPath(new URL('..', import.meta.url)));

const args = process.argv.slice(2);
function argValue(names, fallback) {
  for (let i = 0; i < args.length; i += 1) {
    const [key, inline] = args[i].split('=');
    if (names.includes(key)) {
      if (inline !== undefined) return inline;
      if (i + 1 < args.length && !args[i + 1].startsWith('-')) return args[i + 1];
      return fallback;
    }
  }
  return fallback;
}

const host = argValue(['--host', '-H'], '127.0.0.1');
const port = Number(argValue(['--port', '-p'], '7100'));

const MIME = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.mjs': 'text/javascript; charset=utf-8',
  '.css': 'text/css; charset=utf-8',
  '.json': 'application/json; charset=utf-8',
  '.txt': 'text/plain; charset=utf-8',
  '.md': 'text/markdown; charset=utf-8',
  '.svg': 'image/svg+xml',
  '.png': 'image/png',
  '.jpg': 'image/jpeg',
  '.ico': 'image/x-icon',
};

const server = createServer(async (request, response) => {
  try {
    const url = new URL(request.url ?? '/', `http://${host}`);
    let pathname = decodeURIComponent(url.pathname);
    if (pathname.endsWith('/')) pathname += 'index.html';
    const filePath = normalize(join(root, pathname));
    if (!filePath.startsWith(root)) {
      response.writeHead(403).end('Forbidden');
      return;
    }
    const info = await stat(filePath).catch(() => null);
    const target = info?.isDirectory() ? join(filePath, 'index.html') : filePath;
    const body = await readFile(target);
    response.writeHead(200, {
      'Content-Type': MIME[extname(target).toLowerCase()] ?? 'application/octet-stream',
      'Cache-Control': 'no-store',
    });
    response.end(body);
  } catch {
    response.writeHead(404, { 'Content-Type': 'text/plain; charset=utf-8' }).end('Not Found');
  }
});

server.listen(port, host, () => {
  console.log(`物理模拟实验台 → http://${host}:${port}`);
});

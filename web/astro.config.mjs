import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';
import react from '@astrojs/react';

export default defineConfig({
  site: 'https://vcfkit.dev',
  vite: {
    build: {
      // /wasm/vcfkit_core.js is served as a static public asset; don't try to bundle it.
      rollupOptions: {
        external: ['/wasm/vcfkit_core.js'],
      },
    },
  },
  integrations: [
    react(),
    starlight({
      title: 'vcfkit',
      description: 'Fast VCF operations — normalize, liftover, filter — as a single static binary with zero dependencies.',
      social: [
        { icon: 'github', label: 'GitHub', href: 'https://github.com/robertlangdonn/vcfkit' },
      ],
      sidebar: [
        { label: 'Introduction', link: '/' },
        { label: 'Install', link: '/install' },
        {
          label: 'Commands',
          items: [
            { label: 'normalize', link: '/commands/normalize' },
            { label: 'liftover', link: '/commands/liftover' },
            { label: 'filter', link: '/commands/filter' },
          ],
        },
        { label: 'Benchmarks', link: '/benchmarks' },
        { label: 'Known differences', link: '/known-differences' },
        { label: 'Privacy', link: '/privacy' },
        { label: 'Credits', link: '/credits' },
      ],
      customCss: ['./src/styles/custom.css'],
      components: {
        Hero: './src/components/Hero.astro',
      },
      head: [
        {
          tag: 'meta',
          attrs: { property: 'og:image', content: 'https://vcfkit.dev/og-image.png' },
        },
      ],
    }),
  ],
});

import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';
import react from '@astrojs/react';

export default defineConfig({
  site: 'https://vcfkit.dev',
  vite: {
    optimizeDeps: {
      include: ['react', 'react-dom', 'react-dom/client'],
    },
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
        { label: 'Introduction', link: '/introduction' },
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
        Hero:        './src/components/Hero.astro',
        SiteTitle:   './src/components/SiteTitle.astro',
        SocialIcons: './src/components/SocialIcons.astro',
        Footer:      './src/components/SiteFooter.astro',
      },
      head: [
        {
          tag: 'link',
          attrs: { rel: 'icon', href: '/favicon.svg', type: 'image/svg+xml' },
        },
        {
          tag: 'link',
          attrs: { rel: 'icon', href: '/favicon.ico', sizes: '48x48 32x32 16x16' },
        },
        {
          tag: 'link',
          attrs: { rel: 'apple-touch-icon', href: '/apple-touch-icon.png' },
        },
        {
          tag: 'meta',
          attrs: { property: 'og:image', content: 'https://vcfkit.dev/og-image.png' },
        },
      ],
    }),
  ],
});
